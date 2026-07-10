#!/usr/bin/env python3
"""
V1 gate: finite-difference verification of thermalTopO adjoint sensitivities.

Protocol (docs/derivation.md section 5):
 1. Full base run (primal + adjoint + sensitivities) -> J0, topologySensas1.
 2. For each sampled design cell: central differences on beta_i with the
    primal-only configuration -> FD_i = (J+ - J-)/(2 eps).
 3. Compare FD_i against sens_i; report per-cell ratio and summary stats.
"""

import os
import re
import shutil
import subprocess
import random

CASE = os.path.dirname(os.path.abspath(__file__))
EPS = 5.0e-3
NSAMPLE_PER_GROUP = 8
random.seed(42)

os.chdir(CASE)


def run_solver(logname):
    with open(logname, "w") as lf:
        return subprocess.run(
            ["adjointOptimisationFoam"], stdout=lf, stderr=subprocess.STDOUT)


def read_scalar_field(path, n_expected=None):
    """Parse an ascii OpenFOAM volScalarField internalField."""
    txt = open(path).read()
    m = re.search(
        r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(', txt)
    if m:
        n = int(m.group(1))
        start = m.end()
        vals = txt[start:].split(')')[0].split()
        assert len(vals) == n, f"{path}: {len(vals)} != {n}"
        return [float(v) for v in vals]
    m = re.search(r'internalField\s+uniform\s+([0-9eE+.\-]+)', txt)
    if m:
        assert n_expected, "uniform field needs n_expected"
        return [float(m.group(1))] * n_expected
    raise RuntimeError(f"cannot parse {path}")


def write_alpha(vals):
    body = "\n".join(f"{v:.12g}" for v in vals)
    content = f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object alpha; }}
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar>
{len(vals)}
(
{body}
)
;
boundaryField
{{
    ".*" {{ type zeroGradient; }}
    frontAndBack {{ type empty; }}
}}
"""
    open("0/alpha", "w").write(content)


def get_J(log, cap=6000):
    vals = re.findall(r'meanT : ([0-9eE+.\-]+)', open(log).read())
    vals = [v for v in vals if float(v) > 1.0]  # drop weighted objective prints
    if not vals:
        raise RuntimeError(f"no meanT in {log}")
    if len(vals) >= cap:
        raise RuntimeError(
            f"{log}: primal hit iteration cap ({len(vals)}) - unconverged")
    if abs(float(vals[-1]) - float(vals[-20])) > 1e-6:
        raise RuntimeError(
            f"{log}: J still moving at exit "
            f"({vals[-20]} -> {vals[-1]}) - unconverged")
    return float(vals[-1])


def clean_run_dirs():
    for d in ("1", "2", "optimisation"):
        shutil.rmtree(d, ignore_errors=True)


def use_dict(which):  # 'full' or 'primalonly'
    shutil.copy(f"system/optimisationDict.{which}", "system/optimisationDict")


# ---------------------------------------------------------------- prep dicts
full = open("system/optimisationDict.full").read()
primalonly = full.replace(
    """            as1
            {
                active                 true;""",
    """            as1
            {
                active                 false;""")
open("system/optimisationDict.primalonly", "w").write(primalonly)


# ------------------------------------------------------------------ base run
print("[base] full primal+adjoint run", flush=True)
use_dict("full")
clean_run_dirs()
run_solver("log.fd_base")
J0 = get_J("log.fd_base")
sens = read_scalar_field("1/topologySensas1")
shutil.copy("1/topologySensas1", "fd_topologySens.base")
print(f"[base] J0 = {J0:.12g}, ncells = {len(sens)}", flush=True)

alpha0 = read_scalar_field("0/alpha", n_expected=len(sens))
cx = read_scalar_field("0/Cx", n_expected=len(sens))
cy = read_scalar_field("0/Cy", n_expected=len(sens))

# design cells: inside the designSpace box, safely interior in [0,1]
design = [i for i in range(len(sens))
          if 0.02 < cx[i] < 0.08 and 0.002 < cy[i] < 0.018]
print(f"[design] {len(design)} perturbable cells", flush=True)

EPS = 0.008

def make_direction(kind, seed=None):
    d = [0.0]*len(sens)
    if kind == "uniform_blob":
        for i in design:
            if 0.03 < cx[i] < 0.06 and 0.004 < cy[i] < 0.016:
                d[i] = 1.0
    elif kind == "uniform_band":
        for i in design:
            if 0.062 < cx[i] < 0.075:
                d[i] = 1.0
    else:
        rnd = random.Random(seed)
        for i in design:
            d[i] = rnd.uniform(-1.0, 1.0)
    return d

directions = [("rand%d" % k, make_direction("rand", seed=100 + k))
              for k in range(5)]
directions += [("uniform_blob", make_direction("uniform_blob")),
               ("uniform_band", make_direction("uniform_band"))]

use_dict("primalonly")
rows = []
for name, d in directions:
    Js = {}
    for pm, sgn in (("+", 1.0), ("-", -1.0)):
        vals = [alpha0[i] + sgn*EPS*d[i] for i in range(len(sens))]
        write_alpha(vals)
        clean_run_dirs()
        run_solver(f"log.fd_{name}{pm}")
        Js[pm] = get_J(f"log.fd_{name}{pm}")
    fd = (Js["+"] - Js["-"]) / (2*EPS)
    pred = sum(sens[i]*d[i] for i in design)
    ratio = pred/fd if fd != 0 else float("nan")
    rows.append((name, Js["+"], Js["-"], fd, pred, ratio))
    print(f"[{name}] J+={Js['+']:.12g} J-={Js['-']:.12g} "
          f"FD={fd:.6g} adjoint={pred:.6g} ratio={ratio:.6g}", flush=True)

write_alpha(alpha0)
use_dict("full")

with open("fd_results.csv", "w") as f:
    f.write("direction,Jplus,Jminus,FD,adjointPred,pred_over_FD\n")
    for r_ in rows:
        f.write(",".join(str(v) for v in r_) + "\n")

# least-squares single scale: minimise sum (pred - c*FD)^2
num = sum(r_[4]*r_[3] for r_ in rows)
den = sum(r_[3]*r_[3] for r_ in rows)
c = num/den if den else float("nan")
print("\n=== SUMMARY (least-squares scale) ===")
print(f"J0 = {J0:.12g}")
print(f"c* = {c:.6g}")
wsum, esum = 0.0, 0.0
for r_ in rows:
    name, fd, pred = r_[0], r_[3], r_[4]
    rel = (pred/c - fd)/fd if fd else float("nan")
    print(f"  {name:14s} FD={fd:10.4g}  pred/c*={pred/c:10.4g}  rel.err={rel:+.2%}")
    wsum += abs(fd); esum += abs(pred/c - fd)
print(f"|FD|-weighted mean |rel.err| = {esum/wsum:.2%}")
print("PASS" if esum/wsum < 0.05 else "INVESTIGATE")
print("DONE")

#!/usr/bin/env python3
"""Single-cell FD on high-signal cells: per-cell adjoint accuracy."""

import os
import re
import shutil
import subprocess
import statistics

CASE = os.path.dirname(os.path.abspath(__file__))
EPS = 0.02
os.chdir(CASE)


def run_solver(logname):
    with open(logname, "w") as lf:
        return subprocess.run(
            ["adjointOptimisationFoam"], stdout=lf, stderr=subprocess.STDOUT)


def read_scalar_field(path, n=None):
    txt = open(path).read()
    m = re.search(
        r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(', txt)
    if m:
        return [float(v) for v in txt[m.end():].split(')')[0].split()]
    m = re.search(r'internalField\s+uniform\s+([0-9eE+.\-]+)', txt)
    return [float(m.group(1))]*n


def write_alpha(vals):
    body = "\n".join(f"{v:.12g}" for v in vals)
    open("0/alpha", "w").write(f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object alpha; }}
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
""")


def get_J(log, cap=6000):
    vals = [v for v in re.findall(r'meanT : ([0-9eE+.\-]+)', open(log).read())
            if float(v) > 1.0]
    if len(vals) >= cap:
        raise RuntimeError(f"{log}: hit cap")
    return float(vals[-1])


def clean():
    for d in ("1", "2", "optimisation"):
        shutil.rmtree(d, ignore_errors=True)


sens = read_scalar_field("fd_topologySens.base")
alpha0 = read_scalar_field("0/alpha", len(sens))
cx = read_scalar_field("0/Cx", len(sens))
cy = read_scalar_field("0/Cy", len(sens))

# rank design cells by |sens|; pick top cells from each regime
blob, edge, fluid = [], [], []
for i in range(len(sens)):
    if not (0.02 < cx[i] < 0.08 and 0.002 < cy[i] < 0.018):
        continue
    a = alpha0[i]
    if a > 0.3:
        (blob if (0.033 < cx[i] < 0.057 and 0.005 < cy[i] < 0.015)
         else edge).append(i)
    elif a > 0.01:
        fluid.append(i)

picks = []
for pool, npick, name in ((blob, 8, "blob"),):
    pool.sort(key=lambda i: -abs(sens[i]))
    picks += [(i, name) for i in pool[:npick]]

shutil.copy("system/optimisationDict.primalonly", "system/optimisationDict")
rows = []
for k, (i, grp) in enumerate(picks):
    Js = {}
    for pm, sgn in (("+", 1.0), ("-", -1.0)):
        vals = list(alpha0)
        vals[i] += sgn*EPS
        write_alpha(vals)
        clean()
        run_solver(f"log.fdc_{i}{pm}")
        Js[pm] = get_J(f"log.fdc_{i}{pm}")
    fd = (Js["+"] - Js["-"])/(2*EPS)
    rows.append((i, grp, sens[i], fd))
    print(f"[{k+1}/{len(picks)}] cell {i} ({grp}) sens={sens[i]:.6g} "
          f"FD={fd:.6g}", flush=True)

write_alpha(alpha0)
shutil.copy("system/optimisationDict.full", "system/optimisationDict")

# least-squares scale over these cells, then per-cell errors
num = sum(r[2]*r[3] for r in rows)
den = sum(r[3]*r[3] for r in rows)
c = num/den
print(f"\n=== per-cell (c* = {c:.6g}) ===")
errs = []
for i, grp, s_, fd in rows:
    rel = (s_/c - fd)/fd if fd else float("nan")
    errs.append(abs(rel))
    print(f"  cell {i:5d} {grp:6s} FD={fd:9.4g} pred/c*={s_/c:9.4g} "
          f"err={rel:+.1%}")
print(f"median |err| = {statistics.median(errs):.1%}, "
      f"max |err| = {max(errs):.1%}")
with open("fd_cells.csv", "w") as f:
    f.write("cell,group,sens,FD\n")
    for r in rows:
        f.write(",".join(str(v) for v in r) + "\n")
print("DONE")

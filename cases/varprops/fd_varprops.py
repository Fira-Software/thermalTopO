#!/usr/bin/env python3
"""FD verification of the adjoint gradient with temperature-dependent fluid
AND solid properties active: rho(T), cp(T), k_f(T), mu(T), D_s(T).

Same harness as cases/fdcheck/fd_kT.py (patch p-norm objective, production
configuration with regularisation), so the numbers are directly comparable
with the constant-property and D_s(T)-only campaigns.

The adjoint uses the frozen-property linearisation: dC/dT, dD/dT and dnu/dT
are omitted from the adjoint operator. This campaign measures what that costs.
"""
import os, re, shutil, subprocess, random, statistics

os.chdir(os.path.dirname(os.path.abspath(__file__)))
EPS = 0.008

def run_solver(log):
    with open(log, "w") as lf:
        return subprocess.run(["adjointOptimisationFoam"],
                              stdout=lf, stderr=subprocess.STDOUT)

def readf(path, n=None):
    txt = open(path).read()
    m = re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(', txt)
    if m:
        return [float(v) for v in txt[m.end():].split(')')[0].split()]
    m = re.search(r'internalField\s+uniform\s+([0-9eE+.\-]+)', txt)
    return [float(m.group(1))]*n

def write_alpha(vals):
    body = "\n".join(f"{v:.12g}" for v in vals)
    open("0/alpha","w").write(f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object alpha; }}
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
    vals = [v for v in re.findall(r'pnormT : ([0-9eE+.\-]+)', open(log).read())
            if float(v) > 1.0]
    if not vals:
        raise RuntimeError(f"no pnormT in {log}")
    if len(vals) >= cap:
        raise RuntimeError(f"{log}: hit cap")
    if abs(float(vals[-1]) - float(vals[-2])) > 1e-6:
        raise RuntimeError(f"{log}: J still moving")
    return float(vals[-1])

def clean():
    for d in ("1","2","optimisation"):
        shutil.rmtree(d, ignore_errors=True)

def use(dct):
    shutil.copy(f"system/optimisationDict.{dct}", "system/optimisationDict")

# base full run (primal + adjoint + sensitivities)
use("varprops"); clean(); run_solver("log.vp_base")
J0 = get_J("log.vp_base")
sens = readf("1/topOSensas1")
shutil.copy("1/topOSensas1", "fd_topologySens.varprops")
alpha0 = readf("0/alpha", len(sens))
cx = readf("0/Cx", len(sens)); cy = readf("0/Cy", len(sens))
print(f"[base] J0 = {J0:.12g}", flush=True)

design = [i for i in range(len(sens)) if 0.02 < cx[i] < 0.08 and 0.002 < cy[i] < 0.018]

def direction(kind, seed=None):
    d = [0.0]*len(sens)
    if kind == "blob":
        for i in design:
            if 0.03 < cx[i] < 0.06 and 0.004 < cy[i] < 0.016: d[i] = 1.0
    elif kind == "band":
        for i in design:
            if 0.062 < cx[i] < 0.075: d[i] = 1.0
    else:
        rnd = random.Random(seed)
        for i in design: d[i] = rnd.uniform(-1.0, 1.0)
    return d

dirs = [("blob", direction("blob")), ("band", direction("band")),
        ("rand_a", direction("r", 11)), ("rand_b", direction("r", 22)),
        ("rand_c", direction("r", 33))]

use("varprops_primalonly")
rows = []
for name, d in dirs:
    Js = {}
    for pm, sg in (("+",1.0),("-",-1.0)):
        write_alpha([alpha0[i] + sg*EPS*d[i] for i in range(len(sens))])
        clean(); run_solver(f"log.vp_{name}{pm}")
        Js[pm] = get_J(f"log.vp_{name}{pm}")
    fd = (Js["+"]-Js["-"])/(2*EPS)
    pred = sum(sens[i]*d[i] for i in design)
    rows.append((name, fd, pred))
    print(f"[{name}] FD={fd:.6g} adjoint={pred:.6g}", flush=True)

# single cells: 8 blob interior, 6 edge, 4 sponge (away from box-boundary rows)
blob, edge, sponge = [], [], []
for i in design:
    a = alpha0[i]
    if a > 0.3:
        (blob if (0.033 < cx[i] < 0.057 and 0.005 < cy[i] < 0.015) else edge).append(i)
    elif a > 0.01 and 0.004 < cy[i] < 0.016 and 0.024 < cx[i] < 0.076:
        sponge.append(i)
for pool in (blob, edge, sponge):
    pool.sort(key=lambda i: -abs(sens[i]))
picks = [(i,"blob") for i in blob[:8]] + [(i,"edge") for i in edge[:6]] \
      + [(i,"sponge") for i in sponge[:4]]

cells = []
for k,(i,grp) in enumerate(picks):
    Js = {}
    for pm, sg in (("+",1.0),("-",-1.0)):
        vals = list(alpha0); vals[i] += sg*EPS
        write_alpha(vals); clean(); run_solver(f"log.vpc_{i}{pm}")
        Js[pm] = get_J(f"log.vpc_{i}{pm}")
    fd = (Js["+"]-Js["-"])/(2*EPS)
    cells.append((i,grp,sens[i],fd))
    print(f"[{k+1}/{len(picks)}] cell {i} ({grp}) sens={sens[i]:.6g} FD={fd:.6g}", flush=True)

write_alpha(alpha0); use("varprops")

# least-squares scale from the interior cells (cleanest regime), then errors
num = sum(c[2]*c[3] for c in cells if c[1]=="blob")
den = sum(c[3]*c[3] for c in cells if c[1]=="blob")
c_blob = num/den
print(f"\n=== single cells (c* fitted on blob cells = {c_blob:.6g}) ===")
errs = {}
signs = {"agree": 0, "disagree": 0}
for i,grp,s_,fd in cells:
    rel = (s_/c_blob - fd)/fd if fd else float("nan")
    errs.setdefault(grp, []).append(abs(rel))
    signs["agree" if (s_*fd) > 0 else "disagree"] += 1
    print(f"  cell {i:5d} {grp:6s} FD={fd:9.4g} pred/c*={s_/c_blob:9.4g} err={rel:+.1%}")
for grp, e in errs.items():
    print(f"{grp}: median |err| = {statistics.median(e):.1%}, max = {max(e):.1%}")
print(f"sign agreement: {signs['agree']}/{sum(signs.values())} cells")
print("\n=== directions (same c*) ===")
for name, fd, pred in rows:
    print(f"  {name:8s} FD={fd:9.4g} pred/c*={pred/c_blob:9.4g} "
          f"err={(pred/c_blob-fd)/fd:+.1%}" if fd else f"  {name}: FD~0")
print("DONE")

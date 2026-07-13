#!/usr/bin/env python3
"""FD step-size sweep: find the eps plateau where the finite difference is
neither swamped by primal-convergence bias (eps too small) nor by projection
nonlinearity (eps too large). Only inside that plateau is FD a valid gate."""
import os,re,shutil,subprocess
os.chdir(os.path.dirname(os.path.abspath(__file__)))
def run(log):
    with open(log,"w") as f: subprocess.run(["adjointOptimisationFoam"],stdout=f,stderr=subprocess.STDOUT)
def readf(p):
    t=open(p).read(); m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(',t)
    return [float(v) for v in t[m.end():].split(')')[0].split()]
def write_alpha(v):
    b="\n".join(f"{x:.12g}" for x in v)
    open("0/alpha","w").write(f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object alpha; }}
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar>
{len(v)}
(
{b}
)
;
boundaryField
{{
    ".*" {{ type zeroGradient; }}
    frontAndBack {{ type empty; }}
}}
""")
def getJ(log):
    v=[float(x) for x in re.findall(r'peakT : ([0-9eE+.-]+)',open(log).read()) if float(x)>1]
    if abs(v[-1]-v[-2])>1e-7: raise RuntimeError("J still moving")
    return v[-1]
def clean():
    for d in ("1","2","optimisation"): shutil.rmtree(d,ignore_errors=True)

shutil.copy("system/optimisationDict.sens.B","system/optimisationDict")
clean(); run("log.sw_sens")
sens=readf("1/topOSensas1"); alpha0=readf("0/alpha")
shutil.copy("system/optimisationDict.baseline.B","system/optimisationDict")

CELLS=[628,626]
print(f"{'cell':>5} {'eps':>7} {'FD dJ/dalpha':>14} {'adj/FD':>12}")
for i in CELLS:
    for eps in (0.001,):
        Js={}
        for tag,sg in (("+",1.),("-",-1.)):
            v=list(alpha0); v[i]+=sg*eps
            write_alpha(v); clean(); run(f"log.sw_{i}_{eps}{tag}")
            Js[tag]=getJ(f"log.sw_{i}_{eps}{tag}")
        fd=(Js["+"]-Js["-"])/(2*eps)
        r=sens[i]/fd if fd else float('nan')
        print(f"{i:5d} {eps:7.4f} {fd:14.5e} {r:12.4e}",flush=True)
    print()
write_alpha(alpha0)
print("DONE")

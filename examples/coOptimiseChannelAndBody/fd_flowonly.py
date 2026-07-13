#!/usr/bin/env python3
"""FD gate on STOCK upstream OpenFOAM: simple + adjointSimple + PtLosses.
No thermalTopO in the sensitivity path at all. Filter + projection ON.
adjoint/FD must equal the objective weight (1e-6)."""
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
    # read cycle-1 J from the framework's own objective history (unambiguous)
    import glob
    f=glob.glob("optimisation/objective/*/lossesas1")[0]
    rows=[l.split() for l in open(f).read().splitlines() if l and not l.startswith("#")]
    return float(rows[0][1])
def clean():
    for d in ("1","2","optimisation"): shutil.rmtree(d,ignore_errors=True)

shutil.copy("system/optimisationDict.flowonly","system/optimisationDict")
clean(); run("log.fo_sens")
sens=readf("1/topOSensas1"); alpha0=readf("0/alpha")
shutil.copy("system/optimisationDict.flowonly_primal","system/optimisationDict")
print(f"{'cell':>5} {'eps':>7} {'FD dJ/dalpha':>14} {'adj/FD':>12}  (must be 1e-6)")
for i in (628,626):
    for eps in (5e-4,1e-3,2e-3):
        Js={}
        for tag,sg in (("+",1.),("-",-1.)):
            v=list(alpha0); v[i]+=sg*eps
            write_alpha(v); clean(); run(f"log.fo_{i}{tag}")
            Js[tag]=getJ(f"log.fo_{i}{tag}")
        fd=(Js["+"]-Js["-"])/(2*eps)
        r=sens[i]/fd if fd else float('nan')
        print(f"{i:5d} {eps:7.4f} {fd:14.5e} {r:12.4e}",flush=True)
    print()
write_alpha(alpha0)
print("DONE")

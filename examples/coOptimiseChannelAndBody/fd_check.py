#!/usr/bin/env python3
"""Is the adjoint gradient right in THIS configuration?
Perturb single design cells, re-solve the primal, compare central FD against
topOSensas1. Body cells (alpha=0.7) only, so alpha+-eps stays inside [0,1]."""
import os,re,shutil,subprocess
os.chdir(os.path.dirname(os.path.abspath(__file__)))
EPS=0.005
def run(log):
    with open(log,"w") as f:
        subprocess.run(["adjointOptimisationFoam"],stdout=f,stderr=subprocess.STDOUT)
def readf(p):
    t=open(p).read()
    m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(',t)
    return [float(v) for v in t[m.end():].split(')')[0].split()]
def write_alpha(v):
    body="\n".join(f"{x:.12g}" for x in v)
    open("0/alpha","w").write(
f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object alpha; }}
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<scalar>
{len(v)}
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
def getJ(log):
    v=[float(x) for x in re.findall(r'peakT : ([0-9eE+.-]+)',open(log).read()) if float(x)>1]
    return v[-1]
def clean():
    for d in ("1","2","optimisation"): shutil.rmtree(d,ignore_errors=True)

# 1) one adjoint cycle -> sensitivities.
# J is read from the framework's objective history (cycle 1 = pre-update),
# NOT from the log tail, which would pick up the post-update value.
shutil.copy("system/optimisationDict.sens","system/optimisationDict")
clean(); run("log.fdsens")
sens=readf("1/topOSensas1")
alpha0=readf("0/alpha")
cx=readf("0/Cx"); cy=readf("0/Cy")
J0=float(open("optimisation/objective/0/peakTas1").read().splitlines()[1].split()[1])
print(f"[base] J0={J0:.6f}  max|sens|={max(abs(s) for s in sens):.3e}",flush=True)

# body cells only (alpha ~0.7), inside designSpace, largest |sens|
cand=[i for i in range(len(sens))
      if 0.012<cx[i]<0.108 and 0.004<cy[i]<0.020 and 0.5<alpha0[i]<0.9]
cand.sort(key=lambda i:-abs(sens[i]))
picks=cand[:6]

shutil.copy("system/optimisationDict.baseline","system/optimisationDict")
print(f"{'cell':>6} {'alpha':>6} {'adjoint sens':>14} {'FD dJ/dalpha':>14} {'ratio':>10} sign")
for i in picks:
    Js={}
    for tag,sg in (("+",1.),("-",-1.)):
        v=list(alpha0); v[i]+=sg*EPS
        write_alpha(v); clean(); run(f"log.fd_{i}{tag}")
        Js[tag]=getJ(f"log.fd_{i}{tag}")
    fd=(Js["+"]-Js["-"])/(2*EPS)
    r=sens[i]/fd if fd else float('nan')
    ok="SAME" if sens[i]*fd>0 else "** OPPOSITE **"
    print(f"{i:6d} {alpha0[i]:6.2f} {sens[i]:14.4e} {fd:14.4e} {r:10.3e} {ok}",flush=True)
write_alpha(alpha0)
print("DONE")

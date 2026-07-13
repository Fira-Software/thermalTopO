#!/usr/bin/env python3
"""Which cells actually drive topOSens at the probe cells?

topOSens_k = V * sum_i gbeta_i * (dbeta_i/dalpha_k)

We measure dbeta_i/dalpha_k numerically (perturb raw alpha_k, read the solver's
own beta), take gbeta = topologySens (the pre-chain field), and rank the
contributions. Then we look at WHERE the dominant contributors sit: their beta
(grey?), their projection derivative, and crucially whether any of them lie
inside a FIXED zone (which would mean the forward map moves beta in cells the
adjoint chain believes are frozen).
"""
import os,re,shutil,subprocess,math
os.chdir(os.path.dirname(os.path.abspath(__file__)))
EPS=1e-3
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
def clean():
    for d in ("1","2","optimisation"): shutil.rmtree(d,ignore_errors=True)

b,eta=20.0,0.5
den=math.tanh(b*eta)+math.tanh(b*(1-eta))
Pp=lambda x: b*(1-math.tanh(b*(x-eta))**2)/den

shutil.copy("system/optimisationDict.sens","system/optimisationDict")
clean(); run("log.loc_base")
gbeta=readf("1/topologySensas1"); topo=readf("1/topOSensas1")
alpha0=readf("0/alpha"); cx=readf("0/Cx"); cy=readf("0/Cy")

def zone(i):
    if cy[i]>=0.020: return "CRUST(fixed)"
    if cy[i]<=0.004: return "BASE(fixed)"
    if cx[i]<=0.012: return "inPlug(fixed)"
    if cx[i]>=0.108: return "outPlug(fixed)"
    return "design"

shutil.copy("system/optimisationDict.baseline","system/optimisationDict")
for k in (626,628):
    B={}
    for tag,sg in (("+",1.),("-",-1.)):
        v=list(alpha0); v[k]+=sg*EPS
        write_alpha(v); clean(); run(f"log.loc_{k}{tag}")
        B[tag]=readf("1/beta")
    at=readf("1/alphaTilda")
    dbeta=[(p-m)/(2*EPS) for p,m in zip(B["+"],B["-"])]
    contrib=[(gbeta[i]*dbeta[i], i) for i in range(len(dbeta))]
    total=sum(c for c,_ in contrib)
    contrib.sort(key=lambda t:-abs(t[0]))
    print(f"\n=== probe cell {k}  (topOSens={topo[k]:.4e}, V*sum={1.5e-9*total:.4e}) ===")
    print(f"{'cell':>6} {'zone':>13} {'beta':>7} {'P(aT)':>7} {'dbeta/da':>10} {'gbeta':>11} {'contrib':>11} {'%':>6}")
    for c,i in contrib[:8]:
        print(f"{i:6d} {zone(i):>13} {B['+'][i]:7.4f} {Pp(at[i]):7.3f} {dbeta[i]:10.4f} {gbeta[i]:11.3e} {c:11.3e} {100*c/total:6.1f}")
    fixed=sum(c for c,i in contrib if zone(i)!="design")
    print(f"  -> contribution from FIXED zones: {fixed:.3e}  ({100*fixed/total:+.2f}% of total)")
    moved=sum(1 for i in range(len(dbeta)) if abs(dbeta[i])>1e-9 and zone(i)!="design")
    print(f"  -> beta MOVES in {moved} fixed-zone cells under a raw-alpha perturbation")
write_alpha(alpha0)
print("\nDONE")

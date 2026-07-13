#!/usr/bin/env python3
"""CHAIN-ONLY finite difference: no PDE physics involved.

Freeze the pre-chain physical sensitivity  g_beta = dJ/dbeta  (= topologySens).
Define F(alpha) = sum_i g_beta_i * beta_i(alpha), where beta(alpha) is produced
by the solver's OWN filter+projection+fixed-zone machinery.

Then  dF/dalpha_k  must equal the chained sensitivity topOSens_k.

If this fails, the bug is in the chain (projection derivative, filter transpose,
fixed-zone order, volume convention) and NOT in the thermal sensitivity.
"""
import os,re,shutil,subprocess
os.chdir(os.path.dirname(os.path.abspath(__file__)))
EPS=0.001
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

# base run WITH adjoint -> pre-chain (topologySens) and post-chain (topOSens)
shutil.copy("system/optimisationDict.sens","system/optimisationDict")
clean(); run("log.ch_base")
gbeta  = readf("1/topologySensas1")   # dJ/dbeta   (pre-chain)
chained= readf("1/topOSensas1")       # dJ/dalpha  (post-chain)
alpha0 = readf("0/alpha")
beta0  = readf("1/beta")
print(f"max|gbeta(pre-chain)| = {max(abs(x) for x in gbeta):.4e}")
print(f"max|topOSens(chained)|= {max(abs(x) for x in chained):.4e}")
print()

# primal only is enough to regenerate beta(alpha)
shutil.copy("system/optimisationDict.baseline","system/optimisationDict")
print(f"{'cell':>5} {'FD_chain (dF/dalpha)':>21} {'topOSens':>14} {'topOSens/FD_chain':>18} {'beta cells moved':>17}")
for k in (628,626):
    B={}
    for tag,sg in (("+",1.),("-",-1.)):
        v=list(alpha0); v[k]+=sg*EPS
        write_alpha(v); clean(); run(f"log.ch_{k}{tag}")
        B[tag]=readf("1/beta")
    dbeta=[(p-m)/(2*EPS) for p,m in zip(B["+"],B["-"])]
    moved=sum(1 for d in dbeta if abs(d)>1e-9)
    fd_chain=sum(g*d for g,d in zip(gbeta,dbeta))
    r=chained[k]/fd_chain if fd_chain else float('nan')
    print(f"{k:5d} {fd_chain:21.6e} {chained[k]:14.6e} {r:18.4f} {moved:17d}",flush=True)
write_alpha(alpha0)
print("\nEXPECT topOSens/FD_chain == 1.0 if the chain is self-consistent.")

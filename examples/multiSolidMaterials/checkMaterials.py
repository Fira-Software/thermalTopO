#!/usr/bin/env python3
"""Regression: does each labelled zone actually use its own solid diffusivity?"""
import re, sys
def readf(p):
    t=open(p).read()
    m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(',t)
    if m: return [float(v) for v in t[m.end():].split(')')[0].split()]
    u=re.search(r'internalField\s+uniform\s+([0-9eE+.\-]+)',t)
    return None
cx=readf("0/Cx"); cy=readf("0/Cy"); T=readf("1/T"); Ds=readf("1/DSolid")
def zone(i):
    if cx[i]<0.02 and cy[i]<0.004:  return "CuCrZr"
    if cx[i]>0.08 and cy[i]<0.004:  return "copper"
    if 0.02<cx[i]<0.08 and 0.002<cy[i]<0.018: return "designSpace(->tungsten)"
    return "tungsten(default)"
exp={"CuCrZr":4.0e-4, "copper":2.0e-3}
groups={}
for i in range(len(Ds)):
    groups.setdefault(zone(i),[]).append((Ds[i],T[i]))
print(f"{'zone':>26} {'cells':>6} {'mean Ds':>11} {'expected':>11}  result")
ok=True
for z,v in sorted(groups.items()):
    mean=sum(x for x,_ in v)/len(v)
    if z in exp:
        e=exp[z]; good=abs(mean-e)/e < 1e-9
        print(f"{z:>26} {len(v):6d} {mean:11.4e} {e:11.4e}  {'OK' if good else 'MISMATCH'}")
        ok &= good
    else:
        # tungsten is tabulated: D_s(T) = 1.2e-3 - 1e-5*(T-300)
        pred=sum(1.2e-3 - 1e-5*(t-300) for _,t in v)/len(v)
        good=abs(mean-pred)/pred < 1e-9
        print(f"{z:>26} {len(v):6d} {mean:11.4e} {pred:11.4e}  {'OK (tabulated D_s(T))' if good else 'MISMATCH'}")
        ok &= good
print("\nPASS" if ok else "\nFAIL"); sys.exit(0 if ok else 1)

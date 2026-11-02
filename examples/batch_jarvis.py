#!/usr/bin/env python3
"""批量：读 JARVIS-SuperconDB，对全部条目算 T_c，与库值比。"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from zero_supercon import load_jarvis, from_jarvis
from zero_supercon.chain import chain_batch

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "jarvis_supercon.json")

recs = load_jarvis(path)
keep = [r for r in recs
        if r.get("Tc", 0) > 0.05 and r.get("lamb", 0) > 0.02 and (r.get("atoms") or {}).get("elements")]
inputs = from_jarvis(keep, require_tc=False)

# 逐条对齐：from_jarvis 可能跳过个别条目，改走同一条路径逐条算
def one(r):
    from zero_supercon import ChainInput
    from zero_supercon.sources import ATOMIC_MASS
    from zero_supercon.chain import chain, moments_from_a2f
    at = r["atoms"]; els = at.get("elements") or []
    if not els or any(e not in ATOMIC_MASS for e in els):
        return None
    a = float(np.mean(at["abc"]))
    m = float(np.mean([ATOMIC_MASS[e] for e in els]))
    return chain(ChainInput(a=a, M_amu=m, lam=float(r["lamb"]), wlog_K=float(r["wlog"]))).Tc_K

mine, db = [], []
for r in keep:
    t = one(r)
    if t is None or not np.isfinite(t):
        continue
    mine.append(t); db.append(r["Tc"])
mine, db = np.array(mine), np.array(db)
m = (mine > 0) & (db > 0)
ratio = mine[m] / db[m]
print(f"n={m.sum()}  median={np.median(ratio):.4f}  "
      f"[0.5,2]={np.mean((ratio>0.5)&(ratio<2)):.3f}  "
      f"corr(ln)={np.corrcoef(np.log(mine[m]),np.log(db[m]))[0,1]:.4f}")

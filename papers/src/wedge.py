#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""wedge.py —— Zd008：楔形的斜率分解

几何关系（固定 V, M, a）：omega ∝ sqrt(t), lam ∝ V/t  =>  dln omega/dln lam = -1/2.
数据斜率若为 s，则写成 V ∝ t^alpha：s = 0.5/(alpha-1) => alpha = 1 + 0.5/s.
输出：data/wedge.json
"""
from __future__ import annotations
import csv, json, math, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "wedge.json")


def load(path, need_geom=False):
    xs, ys, meta = [], [], []
    for r in csv.DictReader(open(path)):
        try:
            l = float(r["lam"]); w = float(r["wlog_K"])
        except Exception:
            continue
        if l > 0.02 and w > 1:
            xs.append(math.log(l)); ys.append(math.log(w))
            if need_geom:
                meta.append((float(r["M_amu"]), float(r["a_A"])))
    return np.array(xs), np.array(ys), meta


def main():
    lj, wj, meta = load("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv", True)
    lm, wm, _ = load("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv")
    out = {}
    for name, x, y in (("JARVIS", lj, wj), ("MP", lm, wm)):
        s = float(np.polyfit(x, y, 1)[0])
        out[name] = {"n": len(x), "slope": round(s, 3),
                     "corr": round(float(np.corrcoef(x, y)[0, 1]), 3),
                     "alpha": round(1 + 0.5 / s, 3)}
    # JARVIS：固定 M、a 的窗口内的斜率
    meta = np.array(meta)
    bins = []
    for lo, hi in ((10, 20), (20, 40), (40, 80), (80, 160)):
        m = (meta[:, 0] >= lo) & (meta[:, 0] < hi)
        if m.sum() > 30:
            s = float(np.polyfit(lj[m], wj[m], 1)[0])
            bins.append({"M_lo": lo, "M_hi": hi, "n": int(m.sum()), "slope": round(s, 3)})
    out["JARVIS_M_bins"] = bins
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()

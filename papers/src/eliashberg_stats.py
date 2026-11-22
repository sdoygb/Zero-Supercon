#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eliashberg_stats.py —— Zd022 的统计（过滤版）

读 data/eliashberg_full.json，按 Tc_AD>1 K 且比值有限过滤，给分段的谱形效应统计。
输出：data/eliashberg_stats.json
"""
from __future__ import annotations
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = os.path.join(ROOT, "data", "eliashberg_full.json")
OUT = os.path.join(ROOT, "data", "eliashberg_stats.json")


def st(v):
    v = sorted(v)
    q = lambda p: v[min(len(v) - 1, int(p * len(v)))]
    return {"n": len(v), "median": round(q(0.5), 3), "q25": round(q(0.25), 3),
            "q75": round(q(0.75), 3), "min": round(v[0], 3), "max": round(v[-1], 3)}


def main():
    d = json.load(open(IN))
    rows = d["rows"]
    f = [r for r in rows if r["Tc_AD"] > 1.0 and r["ratio"] and 0.05 < r["ratio"] < 20
         and r["shape_effect"] and np.isfinite(r["shape_effect"])]
    bins = [(0.3, 0.6), (0.6, 1.0), (1.0, 1.5), (1.5, 5.0)]
    out = {"n_total": d["n_total"], "n_used": d["n_used"], "n_filtered": len(f),
           "skipped": d["skipped"],
           "shape_effect": st([r["shape_effect"] for r in f]),
           "ratio": st([r["ratio"] for r in f]),
           "by_bin": {}}
    for a, b in bins:
        s = [r["shape_effect"] for r in f if a <= r["lam"] < b]
        out["by_bin"][f"{a}-{b}"] = st(s) if s else {}
    res = [abs(r["shape_effect"] - 1) for r in f]
    xs = [r["wpeak_over_wlog"] for r in f]
    out["corr_absdev_wpeak_over_wlog"] = round(float(np.corrcoef(xs, res)[0, 1]), 2)
    # Debye 基线（同 λ, ω_log；Zd022 §1）
    import sys
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import eliashberg_full as E
    out["debye_base"] = {str(l): round(E.base_of(l), 3) for l in (0.5, 1.0, 2.0)}
    out["bcs_check"] = E.bcs_limit()
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()

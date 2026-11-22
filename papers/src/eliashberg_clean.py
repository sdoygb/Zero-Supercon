#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eliashberg_clean.py —— Zd022 的干净样本：稳定且无虚频的谱

背景：1058 条里 369 条"由谱重算的 λ 与库值差 >10%"，其中 89% 含负频（虚声子）、
      80% 的 stability 是 unstable——不是归一不同，是谱本身不干净。
本脚本只取 stability=='stable' 且无负频的谱，重跑全谱 Eliashberg。

输出：data/eliashberg_clean.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

from constants import C_N  # noqa: F401  (确保 constants 可用)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "eliashberg_clean.json")
SRC = "/Users/oygb/Downloads/Zero/_staging/材料库/jarvis_epc_data_figshare_1058.json"

import sys
sys.path.insert(0, os.path.join(ROOT, "src"))
import eliashberg_full as E  # noqa: E402


def st(v):
    v = sorted(v)
    q = lambda p: v[min(len(v) - 1, int(p * len(v)))]
    return {"n": len(v), "median": round(q(0.5), 3), "q25": round(q(0.25), 3),
            "q75": round(q(0.75), 3), "min": round(v[0], 3), "max": round(v[-1], 3)}


def main():
    d = json.load(open(SRC))
    rows, skipped = [], {"unstable": 0, "imaginary": 0, "bad": 0, "mismatch": 0}
    for r in d:
        if r.get("stability") != "stable":
            skipped["unstable"] += 1
            continue
        x = np.asarray(r["a2F_original_x"], float)
        y = np.asarray(r["a2F_original_y"], float)
        if np.any((x < 0) & (y > 0)):
            skipped["imaginary"] += 1
            continue
        ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y >= 0)
        x, y = x[ok], y[ok]
        if len(x) < 10:
            skipped["bad"] += 1
            continue
        lam, wlog = E.moments(x, y)
        if not np.isfinite(lam) or lam <= 0.02 or not np.isfinite(wlog):
            skipped["bad"] += 1
            continue
        if abs(lam - r["lamb"]) / max(r["lamb"], 1e-9) > 0.10:
            skipped["mismatch"] += 1
            continue
        tE = E.tc_full(x, y)
        if not np.isfinite(tE):
            skipped["bad"] += 1
            continue
        tA = E.ad(lam, wlog) * E.MEV_K
        base = E.base_of(lam)
        rows.append({"jid": r["jid"], "lam": round(lam, 4), "wlog_K": round(wlog * E.MEV_K, 1),
                     "Tc_full": round(tE, 3), "Tc_AD": round(tA, 3),
                     "ratio": round(tE / tA, 3) if tA > 0 else None,
                     "shape_effect": round((tE / tA) / base, 3) if (tA > 0 and base > 0) else None,
                     **E.shape(x, y, lam)})
    f = [r for r in rows if r["Tc_AD"] > 1.0 and r["shape_effect"] and 0.05 < r["ratio"] < 20]
    bins = [(0.1, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 1.5), (1.5, 5.0)]
    out = {"n_source": len(d), "n_clean": len(rows), "n_filtered": len(f), "skipped": skipped,
           "shape_effect": st([r["shape_effect"] for r in f]),
           "by_bin": {f"{a}-{b}": (st([r["shape_effect"] for r in f if a <= r["lam"] < b])
                                   if [r for r in f if a <= r["lam"] < b] else {})
                      for a, b in bins},
           "rows": rows}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("clean", len(rows), "filtered", len(f), "skipped", skipped)
    print("shape", out["shape_effect"])
    for k, v in out["by_bin"].items():
        print(k, v)
    print("saved", OUT)


if __name__ == "__main__":
    main()

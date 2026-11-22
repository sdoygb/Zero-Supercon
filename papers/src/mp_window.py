#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mp_window.py —— 第二张表：MP 声子分档表（7678 条）的常温窗口检查

口径：
  omega_log 来自 phonopy（DFT 声子），lam 来自库内分档模型（F，模型值）。
  窗口：lam >= lam_req(omega; mu*=0.1)，lam_req 见 Zd003。
输出：data/mp_window.json
"""
from __future__ import annotations
import csv, json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"
OUT = os.path.join(ROOT, "data", "mp_window.json")


def num(x):
    try:
        return float(x)
    except Exception:
        return None


def lam_req(w, ms=0.1, T=300.0):
    if w <= 1.2 * T * math.exp(1.04 / (1 - 0.62 * ms)):
        return None
    rhs = math.log(w * (1 / 1.2) / T)
    d = rhs * (1 - 0.62 * ms) - 1.04
    if d <= 0:
        return None
    return (1.04 + rhs * ms) / d


def main():
    rs = []
    for r in csv.DictReader(open(SRC)):
        w, l = num(r["wlog_K"]), num(r["lam"])
        if w and l and w > 1 and l > 0.02:
            rs.append({"w": w, "l": l, "f": r["formula"], "id": r["mpid"]})

    wmin = 1.2 * 300 * math.exp(1.04 / (1 - 0.62 * 0.1))
    hi = [r for r in rs if r["w"] > wmin]
    inside = [r for r in hi if lam_req(r["w"]) and r["l"] >= lam_req(r["w"])]
    closest = sorted(hi, key=lambda r: -(r["l"] / (lam_req(r["w"]) or 1e9)))[:8]
    pareto = []
    for x in (0.2, 0.5, 1.0, 2.0, 5.0):
        sub = [r["w"] for r in rs if r["l"] >= x]
        if sub:
            pareto.append({"lam_ge": x, "n": len(sub), "w_max": max(sub)})

    out = {
        "n": len(rs),
        "w_min_for_300K": round(wmin, 1),
        "n_above_wmin": len(hi),
        "max_w": max(r["w"] for r in rs),
        "max_l_above_wmin": max((r["l"] for r in hi), default=None),
        "n_inside": len(inside),
        "closest": [{"formula": r["f"], "w": round(r["w"], 1), "lam": round(r["l"], 3),
                     "lam_req": round(lam_req(r["w"]), 2),
                     "ratio": round(r["l"] / lam_req(r["w"]), 3)} for r in closest],
        "pareto": pareto,
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "closest"}, ensure_ascii=False, indent=1))
    print("closest:")
    for c in out["closest"]:
        print(" ", c)
    print("saved", OUT)


if __name__ == "__main__":
    main()

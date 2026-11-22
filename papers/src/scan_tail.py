#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scan_tail.py —— Zd010：候选扫描（先过 omega_log 天花板，再谈 lambda）

输出：
  - JARVIS：max wlog、wlog>1091 的条数、lam>2/5 尾部的 max wlog
  - lam0 = lam*(w/300)^2 的前 12 名（"刚度校正耦合"的代理）
  - MP：wlog>1091 的条数与 max lam（复用 mp_window 口径）
数据：data/scan_tail.json
"""
from __future__ import annotations
import csv, json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "scan_tail.json")
WREF = 300.0
WMIN = 300 / (1 / 1.2 * math.exp(-1.04 / (1 - 0.62 * 0.1)))  # 1091 K


def load(path):
    out = []
    for r in csv.DictReader(open(path)):
        try:
            l = float(r["lam"]); w = float(r["wlog_K"])
        except Exception:
            continue
        if l > 0 and w > 1:
            out.append({"l": l, "w": w, "f": r["formula"],
                        "e": r.get("elements", ""),
                        "a": float(r["a_A"]) if "a_A" in r and r["a_A"] else 0.0,
                        "M": float(r["M_amu"]) if "M_amu" in r and r["M_amu"] else 0.0})
    return out


def main():
    j = load("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv")
    m = load("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv")
    top = sorted(j, key=lambda r: -r["l"] * (r["w"] / WREF) ** 2)[:12]
    out = {
        "wmin": round(WMIN, 1),
        "JARVIS": {
            "n": len(j), "max_w": max(r["w"] for r in j),
            "n_above_wmin": sum(1 for r in j if r["w"] > WMIN),
            "lam_gt2_max_w": max(r["w"] for r in j if r["l"] > 2),
            "lam_gt5_max_w": max(r["w"] for r in j if r["l"] > 5),
            "n_lam_gt2": sum(1 for r in j if r["l"] > 2),
            "n_lam_gt5": sum(1 for r in j if r["l"] > 5),
        },
        "MP": {
            "n": len(m), "max_w": max(r["w"] for r in m),
            "n_above_wmin": sum(1 for r in m if r["w"] > WMIN),
            "max_lam_above": max(r["l"] for r in m if r["w"] > WMIN),
        },
        "lam0_top": [{"f": r["f"], "e": r["e"], "lam": round(r["l"], 3),
                      "w": round(r["w"], 1),
                      "lam0": round(r["l"] * (r["w"] / WREF) ** 2, 2),
                      "ceiling_K": round(0.275 * r["w"])} for r in top],
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()

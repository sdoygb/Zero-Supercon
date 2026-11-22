#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""window_scan.py —— Zd025 §3：用修正后的窗口复扫库内两张表

窗口的两端：
  下界（AD 内）：ω_log ≥ 300/0.28 ≈ 1071 K（Zd003 定理 1.1）
  上界修正：强耦合段全谱比 AD 高 1.4–1.5 倍（Zd023），故"实际"门槛可低到
            ω_log ≥ 300/(0.28*1.45) ≈ 739 K——这是乐观端。
本脚本对每张表算 Tc_AD 与 Tc_corrected（λ>1.5 时 ×1.45），报：
  过 300 K 的条数、最大的 Tc_corrected、最大的 ω_log、最好的 10 条、稳定性分布。

输出：data/window_scan.json
"""
from __future__ import annotations
import csv
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "window_scan.json")

JARVIS = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"
MP = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"
KB = 8.617333e-5
BOOST = 1.45


def ad(lam, wlog_K, mu):
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 1e-9:
        return 0.0
    return wlog_K / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def jarvis():
    rows = []
    for r in csv.DictReader(open(JARVIS)):
        try:
            lam, w, a, M = float(r["lam"]), float(r["wlog_K"]), float(r["a_A"]), float(r["M_amu"])
            ms = float(r["mustar"])
            tc = float(r["Tc_K"])
        except Exception:
            continue
        rows.append({"id": r["formula"], "lam": lam, "wlog": w, "mustar": ms, "Tc_tab": tc,
                     "stab": "?"})
    return rows


def mp():
    rows = []
    for r in csv.DictReader(open(MP)):
        try:
            lam, w = float(r["lam"]), float(r["wlog_K"])
        except Exception:
            continue
        rows.append({"id": r["formula"], "lam": lam, "wlog": w, "mustar": 0.10,
                     "Tc_tab": None, "stab": "?"})
    return rows


def scan(rows, tag):
    for r in rows:
        r["Tc_AD"] = ad(r["lam"], r["wlog"], r["mustar"])
        r["Tc_corr"] = r["Tc_AD"] * (BOOST if r["lam"] > 1.5 else 1.0)
    ok = [r for r in rows if r["Tc_corr"] >= 300]
    top = sorted(rows, key=lambda r: -r["Tc_corr"])[:10]
    hi_w = max(rows, key=lambda r: r["wlog"])
    return {"n": len(rows), "n_pass_300": len(ok),
            "max_Tc_AD": round(max(r["Tc_AD"] for r in rows), 1),
            "max_Tc_corr": round(max(r["Tc_corr"] for r in rows), 1),
            "max_wlog": round(hi_w["wlog"], 1),
            "max_wlog_formula": hi_w["id"],
            "top10": [{"id": r["id"], "lam": round(r["lam"], 3), "wlog": round(r["wlog"], 1),
                       "Tc_AD": round(r["Tc_AD"], 1), "Tc_corr": round(r["Tc_corr"], 1)}
                      for r in top]}


def main():
    j = scan(jarvis(), "JARVIS")
    m = scan(mp(), "MP")
    # 窗口的下界与乐观端
    w_ad = 300 / 0.28
    w_opt = 300 / (0.28 * BOOST)
    out = {"boost": BOOST, "wlog_req_ad": round(w_ad), "wlog_req_opt": round(w_opt),
           "JARVIS": j, "MP": m}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("窗口：AD 内 ω_log ≥ %.0f K；乐观端（×1.45）≥ %.0f K" % (w_ad, w_opt))
    for tag, d in (("JARVIS", j), ("MP", m)):
        print("\n%s: n=%d  过 300K=%d  max Tc_AD=%.1f  max Tc_corr=%.1f  max ω_log=%.0f (%s)"
              % (tag, d["n"], d["n_pass_300"], d["max_Tc_AD"], d["max_Tc_corr"],
                 d["max_wlog"], d["max_wlog_formula"]))
        for r in d["top10"][:5]:
            print("   %-10s lam=%.3f wlog=%6.1f  Tc_AD=%6.1f  Tc_corr=%6.1f"
                  % (r["id"], r["lam"], r["wlog"], r["Tc_AD"], r["Tc_corr"]))
    print("\nsaved", OUT)


if __name__ == "__main__":
    main()

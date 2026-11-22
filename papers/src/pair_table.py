#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pair_table.py —— Zd011：层间候选对

A 层：MP 表里 omega_log > 1091 K 的材料（先用 Zd003 的窗口反解 lam_req）
      g_req = sqrt(lam_req * omega_A / (2 N_B))
B 层：JARVIS 表里的金属（V/Nb/Ta/MoN...），用 lam0 = lam*(w/300)^2 作耦合代理
输出：data/pair_table.json
"""
from __future__ import annotations
import csv, json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "pair_table.json")
EVK = 8.617333e-5


def lam_req(w, ms=0.12, T=300.0):
    if w <= 300 / 0.278:
        return None
    rhs = math.log(w * (1 / 1.2) / T)
    d = rhs * (1 - 0.62 * ms) - 1.04
    if d <= 0:
        return None
    return (1.04 + rhs * ms) / d


def main():
    # A candidates
    a = []
    for r in csv.DictReader(open("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv")):
        try:
            w = float(r["wlog_K"]); l = float(r["lam"])
        except Exception:
            continue
        if w > 1091 and l > 0.02:
            lr = lam_req(w)
            if lr is None or lr > 50:
                continue
            ev = w * EVK
            a.append({"f": r["formula"], "w": round(w, 1), "lam": round(l, 3),
                      "lam_req": round(lr, 2),
                      "g": {str(n): round(math.sqrt(lr * ev / (2 * n)), 2) for n in (0.5, 1, 2)},
                      "mpid": r["mpid"]})
    a.sort(key=lambda r: r["g"]["1"])

    # B candidates
    b = []
    for r in csv.DictReader(open("/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv")):
        try:
            l = float(r["lam"]); w = float(r["wlog_K"]); tc = float(r["Tc_K"]) if r["Tc_K"] else 0
        except Exception:
            continue
        if l > 0 and w > 1 and r["formula"] in ("V", "Nb", "Ta", "Mo", "MoN", "Ti", "Zr", "Hf", "W"):
            b.append({"f": r["formula"], "lam": round(l, 3), "w": round(w, 1),
                      "lam0": round(l * (w / 300) ** 2, 2), "Tc": round(tc, 2)})
    b.sort(key=lambda r: -r["lam0"])

    out = {"A": a, "B": b}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("A candidates (sorted by g_req at N=1):")
    for r in a:
        print(f"  {r['f']:10s} w={r['w']:7.1f} lam={r['lam']:.3f} lam_req={r['lam_req']:6.2f} "
              f"g(0.5/1/2)={r['g']['0.5']:.2f}/{r['g']['1']:.2f}/{r['g']['2']:.2f}")
    print("B candidates:")
    for r in b:
        print(f"  {r['f']:4s} lam={r['lam']:.3f} w={r['w']:6.1f} lam0={r['lam0']:.2f} Tc={r['Tc']:.2f}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

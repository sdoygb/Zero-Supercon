#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""layer_route.py —— Zd005 的设计方程：g_req(omega_A, N_B)

lam_req: 解 Allen–Dynes 窗口 1.04(1+l)/(l-mu*(1+0.62l)) = ln(omega*(1/1.2)/T)
g_req  : sqrt(lam_req * omega_A / (2 N_B))
输出：data/layer_route.json
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "layer_route.json")
EV_PER_K = 8.617333e-5


def lam_req(w_K, ms=0.12, T=300.0):
    rhs = math.log(w_K * (1 / 1.2) / T)
    return (1.04 + rhs * ms) / (rhs * (1 - 0.62 * ms) - 1.04)


def main():
    rows = []
    for w in (1400, 1700, 2000, 2500):
        l = lam_req(w)
        for n in (0.5, 1.0, 2.0):
            g = math.sqrt(l * (w * EV_PER_K) / (2 * n))
            rows.append({"omega_K": w, "omega_eV": round(w * EV_PER_K, 5),
                         "lam_req": round(l, 3), "N_B": n, "g_req_eV": round(g, 3)})
    out = {"mu_star": 0.12, "Tc": 300.0, "rows": rows}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"{'w(K)':>6} {'w(eV)':>7} {'lam_req':>8} {'N_B':>5} {'g_req(eV)':>10}")
    for r in rows:
        print(f"{r['omega_K']:6d} {r['omega_eV']:7.4f} {r['lam_req']:8.3f} {r['N_B']:5.1f} {r['g_req_eV']:10.3f}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

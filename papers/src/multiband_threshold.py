#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""multiband_threshold.py —— P2：多带的 λ12 与 μ*12 阈值（模型级）

口径（Zd019）：
  λ_tot = (λ0-μ*0) + N·net,   net = λ12-μ*12
  λ12 = 2 N_B M_eff g_I²/ω_I   （共享界面模，M_eff 个近简并模相干相加）
  μ*12 = μ12/(1+μ12 ln L),    L = ε_F/ω_I

给定 ω_I、g_I、λ0、μ*0、μ12、ε_F，反解：
  - 达到 300 K 需要 λ_req(ω_I)，故 net_req = (λ_req-(λ0-μ*0))/N；
  - M_eff 的阈值 = (net_req+μ*12) / (2 N_B g_I²/ω_I)。

输入分档：g_I 与 ω_I 是模型/材料量（C）；μ12、ε_F 是 C；N、N_B 是设计输入（I）。
输出：data/multiband_threshold.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "multiband_threshold.json")

KB_EV = 8.617333e-5
G_I = 0.0612          # eV（Zd020 §5.5 的界面值）
OMEGA_I_K = 1500.0
LAM0, MU0 = 0.46, 0.12
EF_EV = 5.0


def lam_req(omega_K, mu=0.12, T=300.0):
    """解 T = w/1.2 exp(-1.04(1+l)/(l-mu(1+0.62l)))；超范围返回 inf。"""
    w = omega_K          # K（T_c 与 ω 都用 K）
    a, b = 1.04, mu
    def f(l):
        den = l - b * (1 + 0.62 * l)
        if den <= 1e-9:
            return w / 1.2 - T      # 指数为正无穷 -> Tc=0，函数为 w/1.2
        ex = -a * (1 + l) / den
        if ex > 50:
            return w / 1.2 * math.exp(50) - T
        return w / 1.2 * math.exp(ex) - T
    grid = [0.2 + 0.01 * i for i in range(0, 6000)]
    for i in range(len(grid) - 1):
        if f(grid[i]) * f(grid[i + 1]) <= 0:
            lo, hi = grid[i], grid[i + 1]
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if f(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
    return float("inf")


def main():
    wI = OMEGA_I_K * KB_EV
    L = EF_EV / wI
    lnL = math.log(L)
    out = {"omega_I_K": OMEGA_I_K, "g_I_eV": G_I, "lam0": LAM0, "mu0": MU0,
           "eps_F_eV": EF_EV, "L": round(L, 2), "lnL": round(lnL, 3),
           "lam_req": {}, "mu12_table": {}, "M_req": {}, "single_band": {}}
    for w in (1000.0, 1200.0, 1500.0, 1916.0, 2200.0):
        out["lam_req"][str(w)] = round(lam_req(w), 3)
    for mu12 in (0.05, 0.10, 0.20):
        out["mu12_table"][str(mu12)] = round(mu12 / (1 + mu12 * lnL), 4)
    lam12_single = 2 * G_I ** 2 / wI          # per unit N_B
    out["lam12_per_NB"] = round(lam12_single, 4)
    for mu12 in (0.05, 0.10, 0.20):
        mu12s = mu12 / (1 + mu12 * lnL)
        for N in (1, 20, 40, 80, 160):
            net_req = (out["lam_req"]["1500.0"] - (LAM0 - MU0)) / N
            for NB in (0.2, 0.4, 1.0):
                key = f"mu12={mu12},N={N},N_B={NB}"
                out["M_req"][key] = round(max(0.0, (net_req + mu12s) / (lam12_single * NB)), 1)
    # 单带路线（N=1）与多带路线（N=40）的对照
    for NB in (0.2, 0.4, 1.0):
        net_req1 = out["lam_req"]["1500.0"] - (LAM0 - MU0)
        out["single_band"][f"N_B={NB}"] = round((net_req1 + 0.0732) / (lam12_single * NB), 1)
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "M_req"}, ensure_ascii=False, indent=1))
    print("\nM_eff 阈值（ω_I=1500K, μ12=0.1）:")
    for N in (1, 20, 40, 80, 160):
        row = [out["M_req"][f"mu12=0.1,N={N},N_B={nb}"] for nb in (0.2, 0.4, 1.0)]
        print(f"  N={N:4d}: N_B=0.2/0.4/1.0 -> {row}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

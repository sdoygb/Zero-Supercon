#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""tc_error_budget.py —— Tc 误差预算与界面模频率

两件事：
1) 界面模频率：Morse 键的曲率 K_I=2 a² D_b，ω_I=sqrt(K_I/μ)/2πc。
   对 Ag–C 给 ω_I（这是界面反向散射声子的物理频率，不是石墨烯的 G 模）。
2) 误差预算：把 (D_I, N_B, ω_I, K, μ*) 的区间做蒙特卡洛，
   给 Tc 的中位与 16–84% 带，并逐项给灵敏度 ∂lnTc/∂lnX。

输出：data/tc_error_budget.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "tc_error_budget.json")

AMU = 1.66054e-27
EV_A2 = 1.602176634e-19 / 1e-20      # eV/Å² -> N/m


def omega_I(D_b_eV, a_A, mu_amu):
    K = 2.0 * D_b_eV * a_A ** 2 * EV_A2
    w = math.sqrt(K / (mu_amu * AMU))
    return w / (2 * math.pi * 2.99792458e10), K      # cm^-1, N/m


def tc_ad(lam, w_K, mu):
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 1e-6:
        return 0.0
    return w_K / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def main():
    print("=" * 80)
    print("1) 界面模频率（Morse 曲率）")
    mC, mAg = 12.0, 107.87
    mu_C_Ag = mC * mAg / (mC + mAg)
    rows = []
    for tag, Db, a in (("Ag–C（弱键）", 2.0, 1.8), ("Ag–C（中）", 2.5, 2.0), ("Ag–C（强键）", 3.0, 2.2)):
        cm, K = omega_I(Db, a, mu_C_Ag)
        rows.append({"bond": tag, "D_b_eV": Db, "a_A": a, "K_N_m": round(K),
                     "omega_cm": round(cm), "omega_K": round(cm * 1.438777)})
        print("  %-14s D=%.1f a=%.1f  K=%.0f N/m  ω=%.0f cm^-1 = %.0f K"
              % (tag, Db, a, K, cm, cm * 1.438777))
    # 参照：金刚石 C–C 光学模
    cm_dia = 1332.0
    print("  参照：金刚石 C–C 光学模 1332 cm^-1 = %.0f K；石墨烯 G 模 1580 cm^-1 = %.0f K"
          % (cm_dia * 1.438777, 1580 * 1.438777))

    print("=" * 80)
    print("2) 误差预算（蒙特卡洛，Ag/GO）")
    rng = np.random.default_rng(20261123)
    n = 40000
    # 输入区间（C/I 分档）
    DI = rng.uniform(3.2, 19.5, n)          # 三种标定的极端（C–C 基准…Harrison η=2.5）
    NB = rng.uniform(0.083, 0.150, n)       # E_F 5–9 eV 的自由电子 s 带
    K = rng.uniform(0.083, 0.250, n)        # slab 标定 ±50%
    MU = rng.uniform(0.10, 0.20, n)
    lam = K * NB * DI ** 2
    out = {}
    for tag, w_lo, w_hi in (("物理界面模（600–1200 K）", 600, 1200),
                            ("乐观（C–C 光学 1336 K）", 1336, 1336),
                            ("最乐观（石墨烯 G 2274 K）", 2274, 2274)):
        W = rng.uniform(w_lo, w_hi, n)
        tc = np.array([tc_ad(l, w, m) for l, w, m in zip(lam, W, MU)])
        out[tag] = {"median": round(float(np.median(tc))),
                    "p16": round(float(np.percentile(tc, 16))),
                    "p84": round(float(np.percentile(tc, 84))),
                    "p5": round(float(np.percentile(tc, 5))),
                    "p95": round(float(np.percentile(tc, 95)))}
        print(f"  {tag:24s} 中位 {out[tag]['median']:4d} K  [{out[tag]['p16']:3d}, {out[tag]['p84']:3d}]（16-84 分位）  [{out[tag]['p5']:3d}, {out[tag]['p95']:3d}]（5-95 分位）")

    print("=" * 80)
    print("3) 灵敏度（在中心点：D_I=8.67, N_B=0.107, ω=1000 K, K=0.166, μ*=0.12）")
    DI0, NB0, K0, MU0, W0 = 8.67, 0.107, 0.166, 0.12, 1000.0
    lam0 = K0 * NB0 * DI0 ** 2
    t0 = tc_ad(lam0, W0, MU0)
    sens = {}
    for name, var, h in (("D_I", DI0, 0.02), ("N_B", NB0, 0.02), ("K", K0, 0.02),
                         ("ω", W0, 0.02), ("μ*", MU0, 0.02)):
        if name == "D_I":
            tp = tc_ad(K0 * NB0 * (DI0 * (1 + h)) ** 2, W0, MU0)
            tm = tc_ad(K0 * NB0 * (DI0 * (1 - h)) ** 2, W0, MU0)
        elif name == "N_B":
            tp = tc_ad(K0 * NB0 * (1 + h) * DI0 ** 2, W0, MU0)
            tm = tc_ad(K0 * NB0 * (1 - h) * DI0 ** 2, W0, MU0)
        elif name == "K":
            tp = tc_ad(K0 * (1 + h) * NB0 * DI0 ** 2, W0, MU0)
            tm = tc_ad(K0 * (1 - h) * NB0 * DI0 ** 2, W0, MU0)
        elif name == "ω":
            tp = tc_ad(lam0, W0 * (1 + h), MU0)
            tm = tc_ad(lam0, W0 * (1 - h), MU0)
        else:
            tp = tc_ad(lam0, W0, MU0 * (1 + h))
            tm = tc_ad(lam0, W0, MU0 * (1 - h))
        d = (math.log(tp) - math.log(tm)) / (2 * math.log(1 + h))
        sens[name] = round(d, 2)
        print("  ∂lnTc/∂ln%-3s = %+6.2f" % (name, d))
    print("  （中心点 Tc = %.0f K）" % t0)

    json.dump({"interface_mode": rows,
               "bands": out,
               "sensitivity": sens,
               "inputs": {"D_I": [3.2, 19.5], "N_B": [0.083, 0.150], "K": [0.083, 0.250],
                          "mu_star": [0.10, 0.20], "omega_phys": [600, 1200]}},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

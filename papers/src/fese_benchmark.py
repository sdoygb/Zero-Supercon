#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fese_benchmark.py —— FeSe/SrTiO3 案例的对照计算

做法：用 Zd029 的 Allen–Dynes 形式**反解** λ——给定观测的 Tc 与界面声子频率，
求所需的 λ。然后与 Ag/GO 的链（λ=K N_B D_I²）对照，看两个机制是否同一通道。

输入（F，文献值）：体相 FeSe Tc=8 K；单层 FeSe/STO Tc=65 K（有 109 K 的报道）；
  STO 的 LO 声子 ~90 meV = 1045 K（界面耦合的声子）。
输出：data/fese_benchmark.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "fese_benchmark.json")


def tc_ad(lam, w, mu):
    den = lam - mu * (1 + 0.62 * lam)
    return w / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def lam_req(tc, w, mu):
    lo, hi = 1e-3, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tc_ad(mid, w, mu) < tc:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    rows = []
    cases = [
        ("体相 FeSe", 8.0, 250.0),
        ("单层 FeSe/STO", 65.0, 1045.0),
        ("单层 FeSe/STO（109 K 报道）", 109.0, 1045.0),
    ]
    for name, tc, w in cases:
        for mu in (0.12, 0.20):
            lr = lam_req(tc, w, mu)
            rows.append({"case": name, "Tc_obs": tc, "omega_K": w, "mu_star": mu,
                         "lambda_req": round(lr, 3)})
            print("%-26s Tc=%5.1f K  ω=%4.0f K  μ*=%.2f  → λ_req=%.3f" % (name, tc, w, mu, lr))

    # Ag/GO 的链
    print("\nAg/GO 的链（Zd029）：")
    DI, K = 8.6726, 0.1660
    for EF in (7.0, 5.5):
        NB = 3.0 / (4 * EF)
        lam = K * NB * DI ** 2
        print("  E_F=%.1f: N_B=%.4f λ=%.3f → Tc_AD=%.1f K（ω=1021 K）"
              % (EF, NB, lam, tc_ad(lam, 1021.0, 0.12)))
    # 若界面耦合是 FeSe/STO 量级（λ≈1.0）
    print("  若 λ=1.0（FeSe/STO 量级）：Tc_AD=%.1f K（ω=1021 K）" % tc_ad(1.0, 1021.0, 0.12))

    # 短程键通道给 Fe–O 什么 λ（说明机制不同）
    print("\n若把 Fe–O 界面当短程键通道（Zd024 的链）：")
    A_B = 0.529177
    rFe, rO = 1.26, 0.66
    R0 = rFe + rO
    kFe = 4.0 / (3 * A_B)
    kO = 4.55 / (2 * A_B)
    kk = 0.5 * (kFe + kO)
    nI = kk * R0 - 1
    tFe = 1.05 * (2.48 / R0) ** (kFe * 2.48 - 1)
    tO = 2.0 * (1.21 / R0) ** (kO * 1.21 - 1)
    tI = math.sqrt(tFe * tO)
    DI_FeO = nI * tI
    import sys
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import dband_dos as D
    import numpy as np
    ctr, N = D.dos(n=40, nb=300, lattice="bcc")
    wband = ctr[-1] - ctr[0]
    ef = D.n_ef(ctr, N, 6.0)
    NB_Fe = float(np.interp(ef, ctr, N)) * (wband / 5.0)
    lam_FeO = K * NB_Fe * DI_FeO ** 2
    print("  R0=%.2f n=%.2f t_geo=%.2f D_I=%.2f eV；N_B(Fe,3d6)=%.2f → λ=%.1f"
          % (R0, nI, tI, DI_FeO, NB_Fe, lam_FeO))
    print("  >>> 远大于反解值（%.1f–%.2f）：FeSe/STO 不是短程键通道" % (rows[2]["lambda_req"], rows[3]["lambda_req"]))

    json.dump({"cases": rows,
               "ag_go": {"D_I": DI, "K": K, "omega_K": 1021.0,
                         "lambda_E7": round(K * (3 / 28) * DI ** 2, 3),
                         "lambda_E55": round(K * (3 / 22) * DI ** 2, 3)},
               "FeO_short_range": {"R0": round(R0, 3), "n_I": round(nI, 3),
                                   "D_I_eV": round(DI_FeO, 2), "N_B_Fe": round(NB_Fe, 2),
                                   "lambda": round(lam_FeO, 1)}},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("\nsaved", OUT)


if __name__ == "__main__":
    main()

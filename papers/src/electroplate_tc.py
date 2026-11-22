#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""electroplate_tc.py —— 家用条件（干冰 195 K）下的可电镀/化学沉积候选

链：C–M 键的 D_I（Zd024）→ fcc/bcc d 带 N_B（Zd026 扩展）→ λ = K N_B D_I²
（K 由 Zd025 的 C–W slab 标定：D_I=8.86、N_B=1 → λ=13.03 ⟹ K=0.166）→
单 Einstein 谱（碳侧光学模）走 Zd022 求解器 → Tc_AD / Tc_full。

输入分档：
  C  金属的 (n_d, W_d, t0, R0)、Slater 的 κ、电镀可行性、磁性
  I  碳侧模频（石墨烯 1580 cm^-1 ≈ 2274 K；金刚石 1336 K）、μ*=0.12、N_B=1 的标定
输出：data/electroplate_tc.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "electroplate_tc.json")
import sys
sys.path.insert(0, os.path.join(ROOT, "src"))
import eliashberg_full as E  # noqa: E402
import dband_dos as D  # noqa: E402

A_B = 0.529177
K_LAM = 0.166          # Zd025 的 slab 标定：λ = K N_B D_I²
MU = 0.12

# 金属: (n_d, W_d[eV], t0[eV], R0[Å], 点阵, 磁性, 电镀/化学沉积)
METALS = {
    "Cu": (10.0, 3.5, 1.0, 2.56, "fcc", "抗磁", "硫酸铜镀液，最容易"),
    "Ni": (8.0, 4.5, 0.9, 2.49, "fcc", "铁磁", "瓦特镀液，容易"),
    "Pd": (10.0, 6.0, 1.1, 2.75, "fcc", "顺磁（近铁磁）", "氯化钯镀液，可做"),
    "Pt": (9.0, 8.0, 1.3, 2.77, "fcc", "顺磁", "氯铂酸镀液，可做但贵"),
    "Ag": (10.0, 6.0, 1.0, 2.89, "fcc", "抗磁", "氰化银镀液，可做"),
    "Au": (10.0, 8.0, 1.2, 2.88, "fcc", "抗磁", "氰化金镀液，可做但贵"),
    "Co": (7.0, 4.5, 1.0, 2.51, "fcc", "铁磁", "硫酸钴镀液，容易"),
    "W": (4.0, 8.0, 1.0, 2.74, "bcc", "顺磁", "不能单独电镀（只能 Ni-W／Co-W 合金）"),
    "Mo": (5.0, 7.0, 1.3, 2.73, "bcc", "顺磁", "不能单独电镀"),
    "Ni80W20": (6.6, 5.2, 0.95, 2.54, "fcc", "弱磁/顺磁", "Ni-W 合金镀液，可做"),
}

# 元素: (Z_eff, n*)；κ=Z_eff/(n* a_B)
SLATER = {"Cu": (4.0, 3.0), "Ni": (4.1, 3.0), "Pd": (5.2, 4.0), "Pt": (5.8, 5.0),
          "Ag": (5.0, 4.0), "Au": (6.2, 5.0), "Co": (4.2, 3.0), "W": (6.15, 5.0),
          "Mo": (4.6, 4.0), "Ni80W20": (4.8, 3.4)}
R_COV = {"Cu": 1.32, "Ni": 1.24, "Pd": 1.39, "Pt": 1.36, "Ag": 1.45, "Au": 1.36,
         "Co": 1.26, "W": 1.30, "Mo": 1.29, "Ni80W20": 1.26}
KAPPA_C = 3.07
T0_C, R0_C = 2.5, 1.54


def kappa(m):
    z, n = SLATER[m]
    return z / (n * A_B)


def pair_DI(m):
    kM = kappa(m)
    kk = 0.5 * (KAPPA_C + kM)
    R = R_COV[m] + 0.77
    n = kk * R - 1.0
    tC = T0_C * (R0_C / R) ** (KAPPA_C * R0_C - 1.0)
    t0, R0 = METALS[m][2], METALS[m][3]
    nM = kM * R0 - 1.0
    tM = t0 * (R0 / R) ** nM
    t = math.sqrt(tC * tM)
    return R, kk, n, t, n * t


_cache = {}


N_S_FLOOR = 1.5 / 14.0     # 自由电子 s 带在 E_F=7 eV 的每自旋 DOS：总 DOS 3n/(2E_F) 的一半，C


def n_b(m):
    meta = METALS[m]
    nd, Wd, lat = meta[0], meta[1], meta[4]
    if lat not in _cache:
        _cache[lat] = D.dos(n=40, nb=300, lattice=lat)
    ctr, N = _cache[lat]
    w = ctr[-1] - ctr[0]
    ef = D.n_ef(ctr, N, nd)
    N0 = float(np.interp(ef, ctr, N))
    return max(N0 * (w / Wd), N_S_FLOOR)


def einstein(w0_meV, lam, npts=800):
    wmax = 2.2 * w0_meV
    grid = np.linspace(1e-2, wmax, npts)
    sig = 0.03 * w0_meV
    a2f = lam * w0_meV / 2.0 * np.exp(-0.5 * ((grid - w0_meV) / sig) ** 2) / (sig * math.sqrt(2 * math.pi))
    return grid, a2f


def main():
    rows = []
    for m, meta in METALS.items():
        R, kk, n, t, DI = pair_DI(m)
        NB = n_b(m)
        lam = K_LAM * NB * DI ** 2
        for w0_name, w0 in (("石墨烯 2274 K", 2274 * 0.0861733), ("金刚石 1336 K", 1336 * 0.0861733)):
            grid, a2f = einstein(w0, lam)
            lam_chk, wlog = E.moments(grid, a2f)
            tA = E.ad(lam_chk, wlog, mu=MU) * E.MEV_K
            tE = E.tc_full(grid, a2f, mu=MU, N=160, tmin=20.0, tmax=1500.0)
            rows.append({"metal": m, "magnetic": meta[5], "plating": meta[6],
                         "R0_A": round(R, 2), "n_I": round(n, 2), "t_I": round(t, 2),
                         "D_I": round(DI, 2), "N_B": round(NB, 3), "lam": round(lam, 2),
                         "mode": w0_name, "wlog_K": round(wlog * E.MEV_K),
                         "Tc_AD": round(tA), "Tc_full": round(tE)})
            print("%-8s %-10s D_I=%5.2f N_B=%.2f lam=%5.2f %-12s AD=%4.0f full=%4.0f"
                  % (m, meta[5][:4], DI, NB, lam, w0_name, tA, tE))
    json.dump({"K_lam": K_LAM, "mu_star": MU, "rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

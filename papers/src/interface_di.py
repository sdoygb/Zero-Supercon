#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_di.py —— 界面形变势 D_I 从 C–W 结构来

链：轨道衰减 → 跳跃 t(R) → 对数导数 n_I → D_I=n_I t_I → slab 的 T_c。

1) 轨道衰减常数（Slater 有效核电荷）：kappa = Z_eff/(n* a_B)。
2) 跳跃的距离依赖：t(R)=A R^{-1} e^{-kappa R}（两轨道尾的乘积），故
   n(R) := -dln t/dln R = kappa R - 1。
3) 标定 A：用已知键（C–C 1.54 Å / 2.5 eV；W–W 2.74 Å / 1.0 eV），
   另给 Harrison 型 t = eta hbar^2/(m d^2) 作第三路。
4) C–W：R0=2.07 Å → t_I、n_I → D_I=n_I t_I。
5) 把 D_I 喂回 slab（Zd018 的构造），算 Tc_AD 与 Tc_full（Zd022 求解器）。

输出：data/interface_di.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_di.json")
import sys
sys.path.insert(0, os.path.join(ROOT, "src"))
import eliashberg_full as E  # noqa: E402
import slab_eliashberg as S  # noqa: E402

A_B = 0.529177          # Å
HB2_ME = 7.619964       # hbar²/m_e in eV·Å²（Harrison 型 t=eta*hbar²/(m d²)）

# (元素, 轨道, Z_eff, n*)，Slater 规则
ORB = {"C": (3.25, 2.0), "W": (6.15, 5.0)}


def kappa(elem):
    z, n = ORB[elem]
    return z / (n * A_B)      # Å^-1


def main():
    kC, kW = kappa("C"), kappa("W")
    R_CC, t_CC = 1.54, 2.5
    R_WW, t_WW = 2.74, 1.0
    R_CW = 2.07
    kCW = 0.5 * (kC + kW)

    # 标定 A（t=A R^{-1} e^{-k R}）：A = t R e^{k R}
    A_CC = t_CC * R_CC * math.exp(kC * R_CC)
    A_WW = t_WW * R_WW * math.exp(kW * R_WW)
    t_CW_CC = A_CC / R_CW * math.exp(-kCW * R_CW)
    t_CW_WW = A_WW / R_CW * math.exp(-kCW * R_CW)
    # Harrison 型：t = eta*hbar^2/(m d^2)，eta 取 1.0/1.84/2.5
    t_CW_H = {str(e): round(e * HB2_ME / R_CW ** 2, 3) for e in (1.0, 1.84, 2.5)}

    n_overlap = kCW * R_CW - 1.0
    n_harrison = (3, 4)     # s-d: 0+2+1；p-d: 1+2+1
    t_lo = min(t_CW_CC, t_CW_WW, min(t_CW_H.values()))
    t_hi = max(t_CW_CC, t_CW_WW, max(t_CW_H.values()))
    n_lo, n_hi = min(n_harrison[0], n_overlap), max(n_harrison[1], n_overlap)
    D_lo, D_hi = n_lo * t_lo, n_hi * t_hi

    # 几何平均标定（异核跳跃的常用估计）：t_CW ≈ sqrt(t_CC(R) t_WW(R))
    n_C_hom = kC * R_CC - 1.0
    n_W_hom = kW * R_WW - 1.0
    t_CC_at = t_CC * (R_CC / R_CW) ** n_C_hom
    t_WW_at = t_WW * (R_WW / R_CW) ** n_W_hom
    t_CW_geo = math.sqrt(t_CC_at * t_WW_at)
    D_I_geo = n_overlap * t_CW_geo

    # slab：D_I 扫描（κ=3，N_B=1）
    slab = []
    for DI in (3.5, 5.0, 7.0, 9.0, 10.0):
        wK, g = S.slab_modes(3.0, DI=DI, NB=1.0)
        grid, a2f = S.a2f_from_modes(wK, g, 1.0)
        lam, wlog = E.moments(grid, a2f)
        tA = E.ad(lam, wlog) * E.MEV_K
        tE = E.tc_full(grid, a2f, mu=0.12, N=200, tmin=50.0, tmax=1000.0)
        slab.append({"D_I": DI, "lam": round(lam, 2), "wlog_K": round(wlog * E.MEV_K, 0),
                     "Tc_AD": round(tA, 0), "Tc_full": round(tE, 0)})

    out = {"R_CW_A": R_CW, "kappa_C": round(kC, 2), "kappa_W": round(kW, 2),
           "kappa_CW": round(kCW, 2), "n_overlap": round(n_overlap, 2),
           "n_harrison": list(n_harrison),
           "t_CW_CC_calib": round(t_CW_CC, 3), "t_CW_WW_calib": round(t_CW_WW, 3),
           "t_CW_harrison": t_CW_H,
           "n_hom_check": {"C-C_overlap": round(n_C_hom, 2), "C-C_empirical": 4,
                           "W-W_overlap": round(n_W_hom, 2), "W-W_empirical": 5},
           "t_CW_geo": round(t_CW_geo, 3), "D_I_geo": round(D_I_geo, 1),
           "t_range": [round(t_lo, 2), round(t_hi, 2)],
           "n_range": [round(n_lo, 2), round(n_hi, 2)],
           "D_I_range": [round(D_lo, 1), round(D_hi, 1)],
           "slab": slab,
           "chain": "kappa -> t(R)=A R^-1 e^-kR -> n=kR-1 -> D_I=n t_I -> slab Tc"}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "slab"}, ensure_ascii=False, indent=1))
    print("slab:")
    for r in slab:
        print("  D_I=%.1f  lam=%.2f  wlog=%.0f  Tc_AD=%.0f  Tc_full=%.0f" % (r["D_I"], r["lam"], r["wlog_K"], r["Tc_AD"], r["Tc_full"]))
    print("saved", OUT)


if __name__ == "__main__":
    main()

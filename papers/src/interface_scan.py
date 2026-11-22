#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_scan.py —— Zd025：C–(W,Ta,Hf,Nb,Mo) 界面键的 D_I 与 slab T_c 扫描

A 侧固定为 C（金刚石／碳氮化物／CN 的界面键都是 C–M）；B 侧五个金属。
D_I=n_I t_I 的口径与 Zd024 同：t(R)=A R^{-1}e^{-kappa R}，异核用几何平均标定。
slab 用 Zd018 的构造，按每个 B 金属的质量泛化（m_B），界面弹簧比 kappa=2/3。

输入分档：
  C  Slater 有效电荷、共价半径、同核参考跳跃 t0 与键长 R0
  I  N_B=1、界面弹簧比 kappa（设计输入）
输出：data/interface_scan.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_scan.json")
import sys
sys.path.insert(0, os.path.join(ROOT, "src"))
import eliashberg_full as E  # noqa: E402
import slab_eliashberg as S  # noqa: E402

A_B = 0.529177
HB2_ME = 7.619964

# 元素: (Z_eff, n*, 共价半径 Å, 同核 t0 eV, 同核 R0 Å, 原子质量)
ELEM = {
    "C": (3.25, 2.0, 0.77, 2.50, 1.54, 12.0),
    "W": (6.15, 5.0, 1.30, 1.00, 2.74, 183.84),
    "Ta": (6.00, 5.0, 1.34, 1.10, 2.86, 180.95),
    "Hf": (5.90, 5.0, 1.44, 1.20, 3.13, 178.49),
    "Nb": (4.50, 4.0, 1.34, 1.30, 2.86, 92.91),
    "Mo": (4.60, 4.0, 1.29, 1.30, 2.73, 95.95),
}


def kappa(e):
    z, n = ELEM[e][0], ELEM[e][1]
    return z / (n * A_B)


def t_at(e, R):
    """同核跳跃按 t=t0 (R0/R)^n 外推，n=kappa R0 -1。"""
    _, _, _, t0, R0, _ = ELEM[e]
    n = kappa(e) * R0 - 1.0
    return t0 * (R0 / R) ** n


def pair(A, B):
    kA, kB = kappa(A), kappa(B)
    kk = 0.5 * (kA + kB)
    R = ELEM[A][2] + ELEM[B][2]
    n = kk * R - 1.0
    tA, tB = t_at(A, R), t_at(B, R)
    t = math.sqrt(tA * tB)
    return {"pair": f"{A}–{B}", "R0_A": round(R, 3), "kappa_A": round(kA, 2),
            "kappa_B": round(kB, 2), "kappa": round(kk, 2), "n_I": round(n, 2),
            "t_A_side": round(tA, 3), "t_B_side": round(tB, 3), "t_I": round(t, 3),
            "D_I": round(n * t, 2)}


def main():
    out = {"pairs": [], "slab": []}
    for B in ("W", "Ta", "Hf", "Nb", "Mo"):
        p = pair("C", B)
        out["pairs"].append(p)
        for kap in (2.0, 3.0):
            wK, g = S.slab_modes(kap, DI=p["D_I"], NB=1.0, m_B=ELEM[B][5])
            grid, a2f = S.a2f_from_modes(wK, g, 1.0)
            lam, wlog = E.moments(grid, a2f)
            tA = E.ad(lam, wlog) * E.MEV_K
            tE = E.tc_full(grid, a2f, mu=0.12, N=200, tmin=50.0, tmax=1500.0)
            out["slab"].append({"pair": p["pair"], "kappa": kap, "D_I": p["D_I"],
                                "lam": round(lam, 2), "wlog_K": round(wlog * E.MEV_K, 0),
                                "Tc_AD": round(tA, 0), "Tc_full": round(tE, 0)})
            print("%-8s kappa=%.0f  D_I=%5.2f  lam=%6.2f  wlog=%4.0fK  AD=%4.0f  full=%4.0f"
                  % (p["pair"], kap, p["D_I"], lam, wlog * E.MEV_K, tA, tE))
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("\n键表:")
    for p in out["pairs"]:
        print("  %-8s R0=%.2f  k=%.2f  n=%.2f  t=%.2f  D_I=%.2f"
              % (p["pair"], p["R0_A"], p["kappa"], p["n_I"], p["t_I"], p["D_I"]))
    print("saved", OUT)


if __name__ == "__main__":
    main()

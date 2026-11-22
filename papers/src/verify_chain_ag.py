#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_chain_ag.py —— Ag/GO 全链条独立复算（第二实现）

逐项：D_I → N_B → K → ω_I → λ → Tc_AD → ξ → 方阻↔厚度 → 四端 R 公式 → 升温功率。
"""
from __future__ import annotations
import math

A_B = 0.529177
EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
E = 1.602176634e-19
HBAR = 1.0546e-34
KB = 8.617333e-5
MU = 0.12


def check(name, got, want, tol):
    ok = abs(got - want) <= tol * max(abs(want), 1e-12)
    print("  %-30s got=%10.4f want=%10.4f  %s" % (name, got, want, "OK" if ok else "FAIL"))


def main():
    print("=" * 76)
    print("1) D_I(C–Ag)")
    rC, rAg = 0.77, 1.45
    R0 = rC + rAg
    kC = 3.25 / (2 * A_B)
    kAg = 5.0 / (4 * A_B)
    kk = (kC + kAg) / 2
    n = kk * R0 - 1
    tC = 2.5 * (1.54 / R0) ** (kC * 1.54 - 1)
    tAg = 1.0 * (2.89 / R0) ** (kAg * 2.89 - 1)
    tg = math.sqrt(tC * tAg)
    DI = n * tg
    print("  R0=%.4f  kC=%.4f  kAg=%.4f  kk=%.4f  n=%.4f" % (R0, kC, kAg, kk, n))
    print("  tC=%.4f  tAg=%.4f  t_geo=%.4f  D_I=%.4f" % (tC, tAg, tg, DI))
    check("D_I (script 8.67)", DI, 8.67, 0.005)

    print("=" * 76)
    print("2) N_B：自由电子 s 带（每自旋）")
    for EF in (5.5, 7.0):
        Nsp = 3.0 / (4 * EF)
        print("  E_F=%.1f eV → 每自旋 N=%.4f /eV" % (EF, Nsp))
    check("N_B at E_F=7 (script 0.107)", 3.0 / (4 * 7.0), 0.107, 0.01)

    print("=" * 76)
    print("3) K 标定（Zd025 的 C–W slab：D_I=8.86、N_B=1、λ=13.03）")
    K = 13.03 / 8.86 ** 2
    check("K (script 0.166)", K, 0.166, 0.01)

    print("=" * 76)
    print("4) 界面模 ω_I（Morse）")
    Db, a, mu = 2.5, 2.0, 12 * 107.87 / (12 + 107.87)
    K_I = 2 * Db * a ** 2 * EV_A2
    om = math.sqrt(K_I / (mu * AMU))
    cm = om / (2 * math.pi * 2.99792458e10)
    print("  K_I=%.1f N/m  ω=%.1f cm^-1 = %.1f K" % (K_I, cm, cm * 1.438777))
    check("ω_I (script 1021 K)", cm * 1.438777, 1021, 0.01)

    print("=" * 76)
    print("5) λ 与 Tc_AD")
    for EF in (5.5, 7.0):
        NB = 3.0 / (4 * EF)
        lam = K * NB * DI ** 2
        w = cm * 1.438777
        den = lam - MU * (1 + 0.62 * lam)
        tc = w / 1.2 * math.exp(-1.04 * (1 + lam) / den)
        print("  E_F=%.1f: λ=%.4f  Tc_AD=%.1f K" % (EF, lam, tc))
    check("λ (script 1.34, E_F=7)", K * 0.107 * DI ** 2, 1.34, 0.02)

    print("=" * 76)
    print("6) 相干长度 ξ(Tc=95 K)")
    D = 1.76 * KB * 95
    for vF in (1.33e6, 4.6e5):
        xi = HBAR * vF / (math.pi * D * E) * 1e9
        print("  v_F=%.2e → ξ=%.1f nm（2ξ=%.0f nm）" % (vF, xi, 2 * xi))

    print("=" * 76)
    print("7) 方阻 ↔ 厚度（ρ_Ag=1.59 µΩ·cm）")
    rho = 1.59e-8
    for t in (5, 10, 15, 25, 50):
        rs_bulk = rho / (t * 1e-9)
        rs_eff = 2.5 * rho / (t * 1e-9)
        print("  t=%2d nm → R_s(体)=%.2f Ω/□，R_s(有效 2.5ρ)=%.2f Ω/□" % (t, rs_bulk, rs_eff))
    print("  >>> 注意：脚本里的“Rs 2–5 Ω/□ ↔ 10–25 nm”是按有效电阻率 2.5ρ 才成立")

    print("=" * 76)
    print("8) 四端公式：R[mΩ]=(Vp−Vn)[µV]/(2I[mA])")
    Vp, Vn, I = 2000.0, -2000.0, 1.0
    print("  Vp=%.0f µV, Vn=%.0f µV, I=%.1f mA → R=%.1f mΩ（=2 Ω）" % (Vp, Vn, I, (Vp - Vn) / (2 * I)))

    print("=" * 76)
    print("9) 升温功率（铜 50 g）")
    Q = 50 * (300 - 77) * 0.385
    print("  ΔQ=%.0f J；1 W→%.0f min，2 W→%.0f min" % (Q, Q / 60, Q / 120))


if __name__ == "__main__":
    main()

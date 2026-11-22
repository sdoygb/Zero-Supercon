#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_chain_cugo.py —— 独立复算 Cu/GO 链条（第二实现，不复用 electroplate_tc 的函数）

逐项核：
  1) D_I(C–Cu)：R0、κ、n_I、t_C、t_Cu、t_geo、D_I —— 手推公式独立写
  2) N_B(Cu)：fcc d 带在 n_d=10 的 N(E_F)，与自由电子 s 带的下限
  3) λ = K N_B D_I²：K 由 slab 标定（C–W，D_I=8.86、N_B=1、λ=13.03）
  4) Tc_AD：AD 公式独立实现；Tc_full：解读器口径（N、tmin/tmax、μ*）
  5) 归一化核验：单 Einstein 谱的矩是否复现 λ 与 ω_log
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A_B = 0.529177
MU = 0.12


def check(name, got, want, tol, note=""):
    ok = abs(got - want) <= tol * max(abs(want), 1e-12)
    print("%-34s got=%12.5f want=%12.5f  %s %s" % (name, got, want, "OK " if ok else "FAIL", note))
    return ok


def main():
    print("=" * 88)
    print("1) D_I 链（手推）")
    R_C, R_Cu = 0.77, 1.32
    R0 = R_C + R_Cu
    kC = 3.25 / (2 * A_B)
    kCu = 4.0 / (3 * A_B)
    kk = 0.5 * (kC + kCu)
    nI = kk * R0 - 1.0
    tC = 2.5 * (1.54 / R0) ** (kC * 1.54 - 1.0)
    tCu = 1.0 * (2.56 / R0) ** (kCu * 2.56 - 1.0)
    tG = math.sqrt(tC * tCu)
    DI = nI * tG
    print("   R0=%.4f  kC=%.4f  kCu=%.4f  kk=%.4f  nI=%.4f" % (R0, kC, kCu, kk, nI))
    print("   tC=%.4f  tCu=%.4f  t_geo=%.4f  D_I=%.4f" % (tC, tCu, tG, DI))
    check("D_I (script 7.53)", DI, 7.53, 0.01)

    print("=" * 88)
    print("2) N_B（fcc d 带 + s 带下限）")
    import sys
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import dband_dos as D
    ctr, N = D.dos(n=40, nb=300, lattice="fcc")
    w = ctr[-1] - ctr[0]
    # 归一化：∫N dE = 5（每自旋 5 轨道）？
    norm = float(np.trapezoid(N, ctr))
    print("   fcc d-DOS 的归一化 ∫N dE = %.4f（应为 5：每自旋 5 轨道）" % norm)
    ef = D.n_ef(ctr, N, 10.0)
    N_d = float(np.interp(ef, ctr, N)) * (w / 9.0)     # Cu 的 W_d=3.5? 这里用 W_d
    print("   d 带在 n_d=10 的 N(E_F)（W_d=3.5 折算）= %.4f /eV" % (float(np.interp(ef, ctr, N)) * (w / 3.5)))
    # 自由电子 s 带：每原子 1 个电子、E_F=7 eV
    #   总 DOS（两自旋）D(E_F) = 3n/(2 E_F)；每自旋 = 一半
    N_s_total = 3.0 / (2.0 * 7.0)
    N_s_spin = N_s_total / 2.0
    print("   s 带：总 DOS=%.4f /eV，每自旋=%.4f /eV" % (N_s_total, N_s_spin))
    print("   >>> 脚本用的 floor = %.3f（是总 DOS 还是每自旋？）" % (3.0 / (2.0 * 7.0)))

    print("=" * 88)
    print("3) λ = K N_B D_I² 与 K 的标定")
    slab = json.load(open(os.path.join(ROOT, "data", "interface_scan.json")))
    row = next(r for r in slab["slab"] if r["pair"] == "C–W" and r["kappa"] == 3.0)
    K = row["lam"] / (1.0 * row["D_I"] ** 2)
    print("   slab 标定：D_I=%.3f、N_B=1 → λ=%.3f ⟹ K=%.4f" % (row["D_I"], row["lam"], K))
    check("K (script 0.166)", K, 0.166, 0.01)
    for name, NB in (("floor=总DOS 0.214", N_s_total), ("floor=每自旋 0.107", N_s_spin)):
        lam = K * NB * DI ** 2
        print("   %-20s λ = %.4f" % (name, lam))

    print("=" * 88)
    print("4) Tc_AD（独立实现）与单 Einstein 谱的矩")
    def tc_ad(lam, wlog, mu=MU):
        den = lam - mu * (1 + 0.62 * lam)
        return wlog / 1.2 * math.exp(-1.04 * (1 + lam) / den)

    for NB, tag in ((N_s_total, "总DOS"), (N_s_spin, "每自旋")):
        lam = K * NB * DI ** 2
        for w0, tagw in ((2274.0, "石墨烯"), (1336.0, "金刚石")):
            print("   %-6s %-6s λ=%.3f  Tc_AD=%.1f K" % (tag, tagw, lam, tc_ad(lam, w0)))
    # 单 Einstein 谱的矩：α²F = λ ω0/2 · δ(ω-ω0) 的 λ 与 ω_log
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import eliashberg_full as E
    for lam_test in (1.0, 2.0):
        w0 = 1336 * 0.0861733
        grid = np.linspace(1e-2, 2.2 * w0, 800)
        sig = 0.03 * w0
        a2f = lam_test * w0 / 2 * np.exp(-0.5 * ((grid - w0) / sig) ** 2) / (sig * math.sqrt(2 * math.pi))
        lm, wl = E.moments(grid, a2f)
        print("   Einstein 谱 λ_in=%.2f → λ_out=%.4f  ω_log=%.1f K" % (lam_test, lm, wl * E.MEV_K))

    print("=" * 88)
    print("5) 自洽核：脚本的 Cu 行与独立值")
    ep = json.load(open(os.path.join(ROOT, "data", "electroplate_tc.json")))
    cu = next(r for r in ep["rows"] if r["metal"] == "Cu" and "金刚石" in r["mode"])
    print("   脚本：D_I=%.2f N_B=%.3f λ=%.2f AD=%.0f full=%.0f" % (cu["D_I"], cu["N_B"], cu["lam"], cu["Tc_AD"], cu["Tc_full"]))
    print("   独立：D_I=%.2f N_B=%.3f λ=%.2f" % (DI, N_s_spin, K * N_s_spin * DI ** 2))


if __name__ == "__main__":
    main()

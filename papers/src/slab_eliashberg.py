#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""slab_eliashberg.py —— 用全谱 Eliashberg 重算 Zd018 的 slab 模型谱

slab（Zd018 §5）：4×4×16 简单立方，8 层 C + 8 层 W，中心弹簧 K_CC=K_WW=1、K_CW=κ；
界面模的耦合 g_i = D_I q u_i，u_i=sqrt(hbar/(2 m_C ω_i)) sqrt(W_i)；λ=2N_B Σ g_i²/ω_i。

本篇要做的事：把同一组模式写成 α²F(ω)=N_B Σ g_i² δ(ω-ω_i)（高斯展宽），
喂给 Zd022 的线性化 Eliashberg 求解器，与 Zd018 的 Allen–Dynes 值比。
slab 的 λ=3–5，正好是 AD 未被验证的强耦合段。

输出：data/slab_eliashberg.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "slab_eliashberg.json")
import sys
sys.path.insert(0, os.path.join(ROOT, "src"))
import eliashberg_full as E  # noqa: E402

KB = 8.617333e-5     # eV/K（能量换算）
KB_J = 1.380649e-23   # J/K（角频率用）
HBAR = 1.0546e-34
AMU = 1.66054e-27
MEV_K = 1.0 / 0.0861733


def slab_modes(kappa, DI=5.0, NB=1.0, wcut=300.0, NC=8, NW=8, Nx=4, Ny=4,
               omega_C=None, m_A=12.0, m_B=184.0):
    """返回 (ω_i[K], g_i[eV], W_i)；与 interface_coupling.slab_tc 同一构造。"""
    Nz = NC + NW
    mCu, mWu, mC = m_A, m_B, m_A * AMU    # 矩阵用 amu（与 Zd018 一致）；u 用 kg
    idx = {}
    atoms = []
    for iz in range(Nz):
        for iy in range(Ny):
            for ix in range(Nx):
                idx[(ix, iy, iz)] = len(atoms)
                atoms.append((ix, iy, iz))
    N = len(atoms)
    mass = np.array([mCu if a[2] < NC else mWu for a in atoms])
    D = np.zeros((3 * N, 3 * N))
    for (ix, iy, iz), i in idx.items():
        for d, (dx, dy, dz) in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            jx, jy, jz = ix + dx, iy + dy, iz + dz
            if dz == 1 and jz >= Nz:
                continue
            jx %= Nx
            jy %= Ny
            j = idx[(jx, jy, jz)]
            K = 1.0 if (iz < NC and jz < NC) else (1.0 if (iz >= NC and jz >= NC) else kappa)
            e = np.zeros(3)
            e[d] = 1.0
            blk = K * np.outer(e, e)
            D[3 * i:3 * i + 3, 3 * i:3 * i + 3] += blk / mass[i]
            D[3 * j:3 * j + 3, 3 * j:3 * j + 3] += blk / mass[j]
            D[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= blk / math.sqrt(mass[i] * mass[j])
            D[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= blk / math.sqrt(mass[i] * mass[j])
    w, v = np.linalg.eigh(D)
    lamC = 4.0 / m_A
    if omega_C is None:
        omega_C = 1336.0 * math.sqrt(12.0 / m_A)
    iface = [idx[(ix, iy, NC - 1)] for ix in range(Nx) for iy in range(Ny)]
    q = 1.0 / 3e-10
    wK, gs = [], []
    for i in range(len(w)):
        wv = omega_C * math.sqrt(max(w[i], 0.0) / lamC)
        if wv <= wcut:
            continue
        W = sum(v[3 * j:3 * j + 3, i] @ v[3 * j:3 * j + 3, i] for j in iface)
        om = wv * KB_J / HBAR
        u = math.sqrt(HBAR / (2 * mC * om)) * math.sqrt(W)
        wK.append(wv)
        gs.append(DI * q * u)
    return np.array(wK), np.array(gs)


def a2f_from_modes(wK, g, NB, npts=1600, sigma_frac=0.04):
    """α²F(ω)=N_B Σ g_i² δ(ω-ω_i)，高斯展宽；全部用 meV（λ=2∫α²F/ω dω 无量纲）。"""
    om = wK * KB * 1e3                      # meV
    g_mev = g * 1e3                         # meV
    wmax = om.max()
    sig = sigma_frac * wmax
    grid = np.linspace(max(1e-2, om.min() - 3 * sig), wmax + 4 * sig, npts)
    a2f = np.zeros_like(grid)
    for o, gi in zip(om, g_mev):
        a2f += NB * (gi ** 2 / 1e3) * np.exp(-0.5 * ((grid - o) / sig) ** 2) / (sig * math.sqrt(2 * math.pi))
    return grid, a2f


def main():
    rows = []
    for kappa in (2.0, 3.0, 5.0, 8.0, 12.0):
        for DI in (5.0, 10.0):
            wK, g = slab_modes(kappa, DI=DI, NB=1.0)
            lam_mode = 2.0 * sum(gi ** 2 / (wi * KB) for wi, gi in zip(wK, g))   # eV
            wlog_mode = math.exp(sum((gi ** 2 / (wi * KB)) * math.log(wi) for wi, gi in zip(wK, g))
                                 / sum(gi ** 2 / (wi * KB) for wi, gi in zip(wK, g)))
            grid, a2f = a2f_from_modes(wK, g, NB=1.0)
            lam, wlog = E.moments(grid, a2f)
            tA = E.ad(lam, wlog) * MEV_K
            tE = E.tc_full(grid, a2f, mu=0.12, N=200, tmin=50.0, tmax=1000.0)
            rows.append({"kappa": kappa, "D_I": DI, "n_modes": int(len(wK)),
                         "lam_modes": round(lam_mode, 3), "lam_a2f": round(lam, 3),
                         "wlog_modes_K": round(wlog_mode, 1), "wlog_a2f_K": round(wlog * MEV_K, 1),
                         "Tc_AD": round(tA, 1), "Tc_full": round(tE, 1),
                         "ratio": round(tE / tA, 3) if tA > 0 else None})
            print("kappa=%.0f DI=%.0f  modes=%d  lam=%.2f (mode %.2f)  wlog=%.0fK (mode %.0f)  "
                  "Tc_AD=%.0f  Tc_full=%.0f  ratio=%.2f"
                  % (kappa, DI, len(wK), lam, lam_mode, wlog * MEV_K, wlog_mode, tA, tE,
                     tE / tA if tA > 0 else float("nan")))
    json.dump({"mu_star": 0.12, "rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dos_constant.py —— 三维半满态密度常数 c_N 的复算（口径统一的第一步）

问题：λ 链与屏蔽链都要 N(0)，但库表与旧 Zd 脚本用的是 0.3713——它不是任何标准归一。
本脚本用两种算法复算，并把采用的常数与它进的两个公式一起落盘。

谱：ε_k=tλ_k-2dt，λ_k=4Σsin²(k_i/2)（Zb005 定理 1.2），t=1、d=3 → ε=-2Σcos k_i，带宽 12。
归一：单自旋态密度 ∫N dε=1。

算法 A（采用）：Monte Carlo。λ_k 的谱测度在 u=k 上均匀（λ=2-2cos u，du/π 恰为均匀测度），
  故 ε=-2Σcos u_i、u_i~U(0,π)。窗口平均 N̄=(1/W)∫_{-W/2}^{W/2}N dε 在 W→0 时给 N(0)。
算法 B（对照）：512³ k 网格 + 直方图。网格法在带心有 +1% 量级的离散偏置（见 W 依赖）。

输出：data/dos_constant.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

from constants import C_LAM_G1, C_LAM_G13, C_N, KAPPA2, N_TOT

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "dos_constant.json")


def mc_dos(W=0.2, M=20_000_000, seed=11):
    """窗口平均的带心 DOS（自旋分辨，归一到 1）。"""
    rng = np.random.default_rng(seed)
    acc = 0
    n = 0
    left = M
    while left:
        m = min(left, 5_000_000)
        u = rng.uniform(0.0, np.pi, (m, 3))
        e = -2.0 * np.cos(u).sum(axis=1)
        acc += int(np.sum(np.abs(e) < 0.5 * W))
        n += m
        left -= m
    return acc / (n * W)


def grid_dos(n=512, nb=2048):
    """k 网格直方图；返回 (N(0), N(0.5))。"""
    c = np.cos(np.arange(n) * 2.0 * np.pi / n)
    lo, hi = -6.0, 6.0
    w = (hi - lo) / nb
    h = np.zeros(nb)
    for i in range(n):
        e = -2.0 * (c[i] + c[:, None] + c[None, :])
        h += np.histogram(e, bins=nb, range=(lo, hi))[0]
    ctr = lo + (np.arange(nb) + 0.5) * w
    dens = h / (n ** 3) / w
    at = lambda x: float(dens[np.argmin(np.abs(ctr - x))])
    return at(0.0), at(0.5)


def main():
    conv = []
    for W in (0.4, 0.2, 0.1, 0.05):
        conv.append({"W": W, "N0": round(mc_dos(W=W, M=8_000_000), 5)})
    n0_mc = conv[1]["N0"]  # W=0.2 的窗口平均
    n0_grid, n05_grid = grid_dos()
    out = {
        "spectrum": "eps_k=-2*sum_i cos k_i (t=1,d=3), band width 12, per spin, int N dE=1",
        "mc_window": conv,
        "n0_mc": round(n0_mc, 5),
        "n0_grid512": round(n0_grid, 5),
        "n0_grid512_at_0.5": round(n05_grid, 5),
        "c_N_adopted": C_N,
        "n_tot": N_TOT,
        "kappa2": round(KAPPA2, 4),
        "C_lam_gL1": round(C_LAM_G1, 4),
        "C_lam_gL13": round(C_LAM_G13, 4),
        "note": "0.3713（旧值）既不是单自旋 0.1427 也不是双自旋 0.2854；库表 mu* 列与旧脚本用它",
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()

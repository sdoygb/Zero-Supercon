#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""constants.py —— Zd 系列共用的口径常数（单一维护点）

口径（2026-11-22 统一；复算见 src/dos_constant.py，账见 data/dos_constant.json）：

  C_N    三维简立方、半满、t=1 的**单自旋**态密度（归一到 1）：N(0)=C_N/t。
         进 λ=N(0)⟨g²⟩/(M⟨ω²⟩)（McMillan 的 N(0) 是单自旋）→ C=C_N·κ²·g_L。
  N_TOT  屏蔽用的**总**态密度 = 2·C_N。
         Thomas–Fermi 的 κ²=4πe²·(dn/dμ) 里 dn/dμ 计两个自旋；
         进 μ=N_TOT·4πê²/(⟨q²⟩+4πê²·N_TOT)（Zb072 命题 11.9 的口径）。

为什么是两个数：同一个 0.3713 曾在两处通用——λ 链要单自旋、屏蔽链要总态密度。
"同一个字指两个对象"是本项目 §3 禁止的，故拆成 C_N 与 N_TOT，屏蔽处一律显式写 2·C_N。

历史：0.3713 是库表（`_staging/批量Tc_全表.csv` 的 μ* 列）与旧 Zd 脚本用的值，
不是任何标准归一（单自旋 0.1427、双自旋 0.2854 都不是它）；账见
`notes/2026-11-18_复盘_推导与计算审计.md` 第 8 项（P1）与
`notes/2026-11-22_口径统一_cN与g.md`。
"""
from __future__ import annotations
import math

# 单自旋态密度（t=1，半满，三维简立方）
# 复算：Monte Carlo 5e9 样（W=0.01）0.142658；512³ 网格直方图 0.14387（+0.9% 离散偏置）
C_N = 0.1427
# 屏蔽用总态密度（两自旋）
N_TOT = 2.0 * C_N
# 费米波矢的格距标度 κ²=(3π²)^{2/3}
KAPPA2 = (3.0 * math.pi ** 2) ** (2.0 / 3.0)

# λ=C·(D/t)² 的两个几何端点（g_L=1 纵声子；g_L=1/3 三维键向平均）
C_LAM_G1 = C_N * KAPPA2          # ≈ 1.366
C_LAM_G13 = C_LAM_G1 / 3.0       # ≈ 0.455

# 单位换算（供 Zd018／Zd020 的 g 与声子频率用）
HBAR = 1.0546e-34                # J·s
AMU = 1.66054e-27                # kg
EV_TO_HZ = 1.0 / 6.582119569e-16  # E[eV] → ν[Hz]

if __name__ == "__main__":
    print(f"C_N     = {C_N}")
    print(f"N_TOT   = {N_TOT}")
    print(f"KAPPA2  = {KAPPA2:.4f}")
    print(f"C(g_L=1)= {C_LAM_G1:.4f}")
    print(f"C(g_L=1/3)={C_LAM_G13:.4f}")

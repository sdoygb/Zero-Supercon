#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_bond_order.py —— 用 Zc036 定义 1.1 对一个真实 C–Ag 界面图算键级 p_ij

图（结构选择 I）：
  石墨烯型 C 团簇（六元环，sp2，每碳 3 个面内邻居）
  + 一个 Ag 原子桥接在一个 C–C 键上（C–Ag 界面键）
Hückel 邻接矩阵：C–C = 1.0，C–Ag = 1.0（同标度，键级是相对量）
占据集 |O| = N/2（Zc036 命题 1.1：半满即电子计数）

输出 p_{ij}(C–Ag)，再按 D=2p|β| → Φ=2Da² → ω（接 Zd029 §4.3）。
"""
from __future__ import annotations
import json, math, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_bond_order.json")

EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
HB = 6.582119569e-16
KB = 8.617333e-5
CM = 1.4387769


def hexagon(n=2):
    """n×n 六元环石墨烯片段（简单做法：n 个环并成一条带）"""
    # 用两个六元环共享一个键（最小界面单元）
    # 顶点：环1(0-5)，环2(4,5,6,7,8,9) 共享 4-5
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
             (4, 6), (6, 7), (7, 8), (8, 9), (9, 5)]
    return edges


def build():
    edges = hexagon()
    nC = 10
    Ag = nC
    # Ag 桥接在 4-5 键上
    edges2 = edges + [(4, Ag), (5, Ag)]
    N = nC + 1
    A = np.zeros((N, N))
    for i, j in edges2:
        A[i, j] = A[j, i] = 1.0
    return A, N, Ag, edges2


def main():
    A, N, Ag, edges = build()
    evals, evecs = np.linalg.eigh(A)
    # Hückel 约定：占据"最高"的 N/2 个本征值（苯核验：键级 +0.667）
    # 碳骨架 10 个 C，每碳 1 个 p 电子；Ag 的 s 电子不计入
    n_e = 10
    occ = n_e // 2
    idx = np.argsort(evals)[-occ:]     # 最高的 occ 个
    print('N=%d occ=%d（碳的 p 电子计数）本征值最小/最大=%.3f/%.3f' % (N, occ, evals[0], evals[-1]))
    P = evecs[:, idx] @ evecs[:, idx].T
    rows = []
    print("图：%d 个 C + 1 个 Ag；边 %d 条；占据 %d/%d" % (N - 1, len(edges), occ, N))
    print("本征值：", np.round(evals, 3))
    for i, j in edges:
        p = 2 * P[i, j]
        is_ag = (i == Ag or j == Ag)
        rows.append({"edge": f"{i}-{j}", "is_interface": is_ag, "p_bond": round(p, 4)})
    print("\n界面键的键级：")
    for r in rows:
        if r["is_interface"]:
            print("  C%d–Ag: p = %.4f" % (int(r["edge"].split("-")[0]), r["p_bond"]))
    pc = [r["p_bond"] for r in rows if not r["is_interface"]]
    print("\nC–C 键级：%s" % np.round(pc, 3))
    p_int = [r["p_bond"] for r in rows if r["is_interface"]]
    pbar = float(np.mean([abs(x) for x in p_int]))
    print("界面键级均值 p = %.4f" % pbar)

    # 接 Zd029 §4.3
    t_I, a, nu = 1.72, 2.33, 3
    mC, mAg = 12.0, 107.87
    D = 2 * pbar * t_I
    K = 2 * D * a * a * EV_A2
    mu = (mC / nu) * mAg / ((mC / nu) + mAg)
    wcm = HB * math.sqrt(K / (mu * AMU)) / KB / CM
    print("\n接 Zd029 §4.3：D=2p|β|=%.2f eV → Φ=%.0f N/m → ω=%.0f cm^-1 (%.0f K)"
          % (D, K, wcm, wcm * CM))
    json.dump({"N": N, "edges": len(edges), "p_interface": round(pbar, 4),
               "rows": rows, "D_eV": round(D, 3), "Phi_Nm": round(K), "omega_cm": round(wcm)},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

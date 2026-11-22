#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""dband_dos.py —— Zd026：bcc 五轨道 d 带的态密度与 N_B

模型：bcc 布拉维格子（原胞基矢 a/2(1,1,-1) 等），最近邻 8 个；
Slater–Koster 的 d-d 双中心积分（ddσ, ddπ, ddδ）给 5×5 跳跃矩阵。
参数取 Andersen 的 canonical 值（ddσ,ddπ,ddδ=(-16.2, 8.75, -1.5)，
只定形状；绝对能量刻度由各金属的 d 带宽 W_d 给（C））。

输出：每个金属的 n_d（C：名义 d 电子数）、N(E_F)/N_avg、N_B=N(E_F)/W_d（eV^-1）。
data/dband_dos.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "dband_dos.json")

DDS, DDP, DDD = -16.2, 8.75, -1.5     # canonical（形状）
# 名义 d 电子数与 d 带宽（C）：(n_d, W_d[eV])
METALS = {"Hf": (2.0, 9.0), "Ta": (3.0, 9.0), "W": (4.0, 8.0),
          "Nb": (4.0, 6.0), "Mo": (5.0, 7.0)}


def sk_dd(l, m, n):
    """方向余弦 (l,m,n) 上的 5×5 d-d 跳跃矩阵（SK 表）。"""
    l2, m2, n2 = l * l, m * m, n * n
    E = np.zeros((5, 5))
    E[0, 0] = 3 * l2 * m2 * DDS + (l2 + m2 - 4 * l2 * m2) * DDP + (l2 * m2 + l2 + m2 - 1) * DDD
    E[1, 1] = 3 * m2 * n2 * DDS + (m2 + n2 - 4 * m2 * n2) * DDP + (m2 * n2 + m2 + n2 - 1) * DDD
    E[2, 2] = 3 * n2 * l2 * DDS + (n2 + l2 - 4 * n2 * l2) * DDP + (n2 * l2 + n2 + l2 - 1) * DDD
    E[3, 3] = 0.75 * (l2 - m2) ** 2 * DDS + 3 * (l2 - m2) ** 2 * DDP + (1 - 0.75 * (l2 - m2) ** 2) * DDD
    E[4, 4] = (n2 - 0.5 * (l2 + m2)) ** 2 * DDS + 3 * n2 * (l2 + m2) * DDP + 0.75 * (l2 + m2) ** 2 * DDD
    E[0, 1] = E[1, 0] = np.sqrt(3) * l * m * m2 * DDS + np.sqrt(3) * l * m * (1 - 2 * m2) * DDP \
        - np.sqrt(3) * l * m * m2 * DDD
    E[0, 2] = E[2, 0] = np.sqrt(3) * l2 * m * n * DDS + np.sqrt(3) * m * n * (l2 - 1) * DDP \
        - np.sqrt(3) * l2 * m * n * DDD
    E[1, 2] = E[2, 1] = np.sqrt(3) * l * m2 * n * DDS + np.sqrt(3) * l * n * (m2 - 1) * DDP \
        - np.sqrt(3) * l * m2 * n * DDD
    E[0, 3] = E[3, 0] = 1.5 * l * m * (l2 - m2) * DDS + 2 * l * m * (m2 - l2) * DDP \
        + 0.5 * l * m * (l2 - m2) * DDD
    E[1, 3] = E[3, 1] = 1.5 * m * n * (l2 - m2) * DDS + m * n * (1 - 2 * (l2 - m2)) * DDP \
        - 0.5 * m * n * (l2 - m2 + 2) * DDD
    E[2, 3] = E[3, 2] = 1.5 * n * l * (l2 - m2) * DDS + n * l * (1 + 2 * (l2 - m2)) * DDP \
        - 0.5 * n * l * (l2 - m2 - 2) * DDD
    E[0, 4] = E[4, 0] = np.sqrt(3) * l * m * (n2 - 0.5 * (l2 + m2)) * DDS \
        + 2 * np.sqrt(3) * l * m * (m2 - l2) * DDP - 0.5 * np.sqrt(3) * l * m * (l2 + m2) * DDD
    E[1, 4] = E[4, 1] = np.sqrt(3) * m * n * (n2 - 0.5 * (l2 + m2)) * DDS \
        + np.sqrt(3) * m * n * (l2 + m2 - 1) * DDP + 0.5 * np.sqrt(3) * m * n * (l2 + m2) * DDD
    E[2, 4] = E[4, 2] = np.sqrt(3) * n * l * (n2 - 0.5 * (l2 + m2)) * DDS \
        + np.sqrt(3) * n * l * (l2 + m2 - 1) * DDP + 0.5 * np.sqrt(3) * n * l * (l2 + m2) * DDD
    return E


def fcc_nn():
    """fcc 最近邻 12 个（单位：a/2）。"""
    dirs = [(1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
            (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
            (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1)]
    out = []
    for sx, sy, sz in dirs:
        v = 0.5 * np.array([sx, sy, sz], float)
        l, m, n = (v / np.linalg.norm(v))
        out.append((v, sk_dd(l, m, n)))
    return out


def bcc_nn():
    """bcc 最近邻 8 个（单位：a/2）。"""
    out = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                v = np.array([sx, sy, sz], float) * 0.5
                nn = np.linalg.norm(v) ** -1
                l, m, n = v / np.linalg.norm(v)
                out.append((v, sk_dd(l, m, n)))
    return out


def dos(n=64, nb=400, lattice="bcc"):
    nn = bcc_nn() if lattice == "bcc" else fcc_nn()
    if lattice == "bcc":
        # bcc 原胞基矢（a=1）：a/2(1,1,-1) 等；倒格矢
        a1 = 0.5 * np.array([1, 1, -1])
        a2 = 0.5 * np.array([1, -1, 1])
        a3 = 0.5 * np.array([-1, 1, 1])
    else:
        a1 = 0.5 * np.array([0, 1, 1])
        a2 = 0.5 * np.array([1, 0, 1])
        a3 = 0.5 * np.array([1, 1, 0])
    A = np.column_stack([a1, a2, a3])
    B = 2 * np.pi * np.linalg.inv(A).T
    frac = np.linspace(0, 1, n, endpoint=False)
    F1, F2, F3 = np.meshgrid(frac, frac, frac, indexing="ij")
    K = (F1[..., None] * B[0] + F2[..., None] * B[1] + F3[..., None] * B[2]).reshape(-1, 3)
    evals = np.zeros((K.shape[0], 5))
    for i, k in enumerate(K):
        H = np.zeros((5, 5), complex)
        for v, T in nn:
            H += T * np.exp(1j * k @ v)
        evals[i] = np.linalg.eigvalsh(H)
    E = evals.ravel()
    hist, edges = np.histogram(E, bins=nb)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    N = hist / (K.shape[0] * (edges[1] - edges[0]))   # 每原子（5 轨道）的态密度
    return ctr, N


def n_ef(ctr, N, nd):
    """给定 d 电子数 nd（0–10），返回费米能（在 5 轨道带内）。"""
    cum = np.cumsum(N) * (ctr[1] - ctr[0])
    cum /= cum[-1]          # 归一 1（10 态）
    return float(np.interp(nd / 10.0, cum, ctr))


def main():
    ctr, N = dos()
    w = ctr[-1] - ctr[0]
    print("band range", round(ctr[0], 2), round(ctr[-1], 2), "width", round(w, 2), "(canonical units)")
    rows = []
    for m, (nd, Wd) in METALS.items():
        ef = n_ef(ctr, N, nd)
        N0 = float(np.interp(ef, ctr, N))
        # 归一化：N 的单位是 1/（canonical 能量）；real 带宽 W_d 映射到 canonical 宽度 w
        # N_real = N_canon * (w / W_d)
        Nb = N0 * (w / Wd)
        rows.append({"metal": m, "n_d": nd, "W_d_eV": Wd,
                     "E_F_norm": round((ef - ctr[0]) / w, 3),
                     "N0_norm": round(N0 * w, 4),
                     "N_B_per_eV": round(Nb, 3)})
        print(m, rows[-1])
    json.dump({"range": [round(ctr[0], 3), round(ctr[-1], 3)], "rows": rows,
               "dos": {"E": [round(x, 3) for x in ctr.tolist()],
                       "N": [round(x, 5) for x in (N * w).tolist()]}},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

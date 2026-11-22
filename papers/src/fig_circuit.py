#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_circuit.py —— 测量电路与 R(T) 示意"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

CU = "#d98a2b"
INS = "#e8e8e8"


def circuit(ax):
    ax.set_xlim(0, 11.5)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # 叠层
    ax.add_patch(Rectangle((4.2, 2.2), 2.2, 3.6, fc=INS, ec="k", lw=1.4))
    ax.text(5.3, 4.0, "叠层\n(N−1 个\n乙醇间隙)", ha="center", va="center", fontsize=9)

    # 电流源：左侧，接顶/底
    ax.add_patch(Circle((1.4, 4.0), 0.55, fc="white", ec="k", lw=1.4))
    ax.text(1.4, 4.0, "I", ha="center", va="center", fontsize=13)
    ax.plot([1.95, 2.7], [4.0, 4.0], "k", lw=1.4)
    ax.plot([2.7, 2.7], [2.2, 5.8], "k", lw=1.4)
    ax.plot([2.7, 4.2], [5.8, 5.8], "k", lw=1.4)
    ax.plot([2.7, 4.2], [2.2, 2.2], "k", lw=1.4)
    ax.annotate("", xy=(3.5, 5.8), xytext=(2.9, 5.8),
                arrowprops=dict(arrowstyle="-|>", color="crimson", lw=1.5))

    # 电压表：右侧，接顶/底（不穿过样品）
    ax.add_patch(Circle((8.4, 4.0), 0.55, fc="white", ec="k", lw=1.4))
    ax.text(8.4, 4.0, "V", ha="center", va="center", fontsize=12)
    ax.plot([6.4, 8.4], [5.6, 5.6], "k", lw=1.4)
    ax.plot([8.4, 8.4], [5.6, 4.55], "k", lw=1.4)
    ax.plot([6.4, 8.4], [2.4, 2.4], "k", lw=1.4)
    ax.plot([8.4, 8.4], [2.4, 3.45], "k", lw=1.4)

    # 恒温器：温度计 + 加热器
    ax.add_patch(Rectangle((9.9, 4.3), 0.7, 0.5, fc="#fddbc7", ec="k"))
    ax.text(10.25, 4.55, "加热", ha="center", va="center", fontsize=7)
    ax.add_patch(Rectangle((9.9, 2.4), 0.7, 1.1, fc="#f4a582", ec="k"))
    ax.text(10.25, 2.95, "Pt100", ha="center", va="center", fontsize=7, rotation=90)
    ax.text(10.25, 1.9, "恒温器", ha="center", fontsize=8)

    # LCR（电容，可选）
    ax.add_patch(Circle((5.3, 1.2), 0.45, fc="white", ec="k", lw=1.0, ls="--"))
    ax.text(5.3, 1.2, "LCR", ha="center", va="center", fontsize=7)

    ax.text(5.75, 7.7, "四端法：电流从顶/底进出，电压从顶/底两侧取", ha="center", fontsize=9)
    ax.text(5.75, 0.5, "R_total = (N−1)·R_gap + R_Cu（R_Cu 可忽略）", ha="center", fontsize=9)


def r_t(ax):
    T = np.linspace(170, 260, 600)
    # 正常态：离子电导随温度下降；200 K 附近转变
    Rn = 3e4 * np.exp(1500 * (1 / T - 1 / 260))
    drop = 1 / (1 + np.exp(-(T - 200) / 1.2))
    R = Rn * drop + 1e-3
    ax.semilogy(T, R, "k-", lw=2)
    ax.axvline(200, color="crimson", ls="--", lw=1)
    ax.text(201, 1e5, "T_BKT", color="crimson", fontsize=9)
    ax.axhline(1e-3, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("R (Ω)")
    ax.set_xlim(170, 260)
    ax.set_ylim(1e-4, 1e7)
    ax.annotate("正常态：离子电导\nR ~ (N−1)·R_gap", xy=(245, 4e4), xytext=(215, 3e5),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.annotate("BKT 以上：\nR ∝ exp(−b/√(T−T_BKT))", xy=(203, 3e2), xytext=(210, 1e3),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.annotate("T_BKT 以下：\nR → 0（I–V 幂律）", xy=(190, 1e-2), xytext=(172, 1e1),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_title("R(T) 示意（非等比）", fontsize=10)


def main():
    fig = plt.figure(figsize=(12, 5.2))
    ax1 = fig.add_axes([0.02, 0.04, 0.56, 0.90])
    ax2 = fig.add_axes([0.63, 0.16, 0.35, 0.72])
    circuit(ax1)
    r_t(ax2)
    fig.suptitle("叠层样品：四端电阻测量与 R(T) 行为", fontsize=13)
    fig.savefig(os.path.join(ROOT, "data", "fig_circuit.png"), dpi=160)
    plt.close(fig)
    print("saved data/fig_circuit.png")


if __name__ == "__main__":
    main()

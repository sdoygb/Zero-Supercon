#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_box.py —— 半导体冷箱（Peltier）里的叠层装置"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

CU = "#d98a2b"
INS = "#e8e8e8"
MEOH = "#bfe3f7"


def main():
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 9)
    ax.axis("off")

    # 冷箱外壳
    ax.add_patch(Rectangle((1.0, 1.0), 8.0, 6.4, fc="#f5f5f5", ec="k", lw=2))
    ax.text(5.0, 8.2, "半导体冷箱（保温层）", ha="center", fontsize=11)

    # Peltier 冷板（底部）＋ 散热器
    ax.add_patch(Rectangle((1.8, 1.35), 6.4, 0.45, fc="#c6dbef", ec="k", lw=1.2))
    ax.text(5.0, 1.58, "Peltier 冷板（TEC）", ha="center", va="center", fontsize=9)
    ax.add_patch(Rectangle((1.8, 0.35), 6.4, 0.6, fc="#a1d99b", ec="k", lw=1.0))
    ax.text(5.0, 0.65, "散热器/风扇", ha="center", va="center", fontsize=8)

    # 叠层（竖立，立在工作台上）
    x0 = 4.2
    for i in range(5):
        x = x0 + i * 0.55
        ax.add_patch(Rectangle((x, 1.8), 0.35, 3.7, fc=INS, ec="k", lw=0.8))
        ax.add_patch(Rectangle((x - 0.05, 1.8), 0.05, 3.7, fc=CU, ec="none"))
        ax.add_patch(Rectangle((x + 0.35, 1.8), 0.05, 3.7, fc=CU, ec="none"))
        if i < 4:
            ax.add_patch(Rectangle((x + 0.42, 1.8), 0.10, 3.7, fc=MEOH, ec="none", alpha=0.85))
    ax.annotate("叠层：N 片，片间乙醇", xy=(6.1, 5.0), xytext=(7.3, 6.9),
                fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))

    # 引脚引出（顶部馈通）
    for i in range(5):
        x = x0 + i * 0.55 + 0.175
        ax.plot([x, x], [5.5, 7.0], color="k", lw=0.8)
    ax.add_patch(Rectangle((4.0, 7.0), 3.2, 0.25, fc="#d9d9d9", ec="k", lw=0.8))
    ax.text(5.6, 7.45, "密封馈通（引脚引出）", ha="center", fontsize=8)

    # Pt100 / 干燥剂
    ax.add_patch(Rectangle((7.7, 4.6), 0.5, 0.8, fc="#f4a582", ec="k", lw=0.8))
    ax.text(7.95, 5.0, "Pt100", ha="center", va="center", fontsize=7, rotation=90)
    ax.add_patch(Rectangle((7.7, 2.2), 0.5, 0.8, fc="#fee0b6", ec="k", lw=0.8))
    ax.text(7.95, 2.6, "干燥剂", ha="center", va="center", fontsize=7, rotation=90)

    # 冷
    ax.annotate("", xy=(2.3, 4.6), xytext=(2.3, 3.0),
                arrowprops=dict(arrowstyle="-|>", color="#2171b5", lw=1.4))
    ax.text(2.0, 3.8, "冷", color="#2171b5", fontsize=9)

    ax.set_title("装置总览：叠层 ＋ 乙醇 ＋ Peltier 冷箱", fontsize=13)
    fig.savefig(os.path.join(ROOT, "data", "fig_box.png"), dpi=160)
    plt.close(fig)
    print("saved data/fig_box.png")


if __name__ == "__main__":
    main()

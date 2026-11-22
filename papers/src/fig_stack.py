#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_stack.py —— 镀铜绝缘片叠层＋乙醇：二维剖面图与三维爆炸图"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Ellipse
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False

CU = "#d98a2b"
INS = "#e8e8e8"
MEOH = "#bfe3f7"


def fig_2d():
    fig = plt.figure(figsize=(11, 7.2))
    ax = fig.add_axes([0.04, 0.06, 0.62, 0.88])
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    plates = [0.9, 2.9, 4.9, 6.9]   # 竖立：沿 x 排开
    y0, y1 = 1.6, 8.4
    ins_w, cu_w = 0.55, 0.10

    for i, x in enumerate(plates):
        # 绝缘片（竖立）
        ax.add_patch(Rectangle((x, y0), ins_w, y1 - y0, fc=INS, ec="k", lw=1.2))
        # 两面铜（左右）
        ax.add_patch(Rectangle((x - cu_w, y0), cu_w, y1 - y0, fc=CU, ec="k", lw=0.6))
        ax.add_patch(Rectangle((x + ins_w, y0), cu_w, y1 - y0, fc=CU, ec="k", lw=0.6))
        # 顶边缘连通
        ax.add_patch(Rectangle((x - cu_w, y1), ins_w + 2 * cu_w, 0.10,
                               fc=CU, ec="k", lw=0.6))
        # 底部触点
        ax.add_patch(Rectangle((x - cu_w, y0 - 0.9), cu_w, 0.9, fc=CU, ec="k", lw=0.6))
        ax.text(x, y0 - 1.1, f"引脚 {i+1}", ha="center", fontsize=8)

    # 乙醇间隙
    for i in range(3):
        xb = plates[i] + ins_w + cu_w
        xt = plates[i + 1] - cu_w
        ax.add_patch(Rectangle((xb, y0), xt - xb, y1 - y0, fc=MEOH, ec="none", alpha=0.85))
        # 绝缘垫片（顶/底封边，兼密封）
        ax.add_patch(Rectangle((xb, y1 - 0.7), xt - xb, 0.7, fc="#9a9a9a", ec="k", lw=0.6))
        ax.add_patch(Rectangle((xb, y0), xt - xb, 0.7, fc="#9a9a9a", ec="k", lw=0.6))
        ax.text((xb + xt) / 2, (y0 + y1) / 2, "乙醇\n(垫片定间隙)", ha="center",
                va="center", fontsize=8, color="#1f5f8b", rotation=90)

    # 电流路径（水平）
    ax.annotate("", xy=(8.6, 4.9), xytext=(0.2, 4.9),
                arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2.2,
                                linestyle="--"))
    ax.text(4.4, 5.15, "I", color="crimson", fontsize=13, fontweight="bold")

    # 标注
    ax.annotate("铜（1–10 μm）", xy=(plates[0] - cu_w / 2, 3.2),
                xytext=(0.2, 2.2), fontsize=9, color=CU,
                arrowprops=dict(arrowstyle="->", color=CU, lw=1))
    ax.annotate("铝板 1 mm（导电、无磁）", xy=(plates[0] + ins_w / 2, 4.2),
                xytext=(0.2, 5.2), fontsize=9, color="#555",
                arrowprops=dict(arrowstyle="->", color="#555", lw=1))
    ax.annotate("铜全包覆（含边缘，两面同电位）", xy=(plates[2] + ins_w / 2, 8.45),
                xytext=(6.2, 9.4), fontsize=9, color=CU,
                arrowprops=dict(arrowstyle="->", color=CU, lw=1))

    # 间距标注
    ax.annotate("", xy=(plates[0] + ins_w + cu_w, 8.9), xytext=(plates[1] - cu_w, 8.9),
                arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
    ax.text((plates[0] + plates[1]) / 2 + 0.2, 9.35, "间距 = 垫片厚",
            ha="center", fontsize=8)

    ax.text(6.55, 1.9, "绝缘垫片（封乙醇）", fontsize=8, color="#555", rotation=90)
    ax.set_title("镀铜铝板叠层：二维剖面（全包覆，每板一引脚，竖立）", fontsize=12)

    # 放大插图：Cu / 乙醇 / Cu
    axi = fig.add_axes([0.70, 0.10, 0.28, 0.34])
    axi.set_xlim(0, 10)
    axi.set_ylim(0, 6)
    axi.axis("off")
    axi.add_patch(Rectangle((0, 0), 10, 1.5, fc=CU, ec="k"))
    axi.add_patch(Rectangle((0, 1.5), 10, 3.0, fc=MEOH, ec="none", alpha=0.85))
    axi.add_patch(Rectangle((0, 4.5), 10, 1.5, fc=CU, ec="k"))
    for x in (2.0, 5.0, 8.0):
        axi.add_patch(Circle((x, 2.2), 0.18, fc="#1f5f8b", ec="none"))
        axi.add_patch(Circle((x, 3.0), 0.18, fc="#1f5f8b", ec="none"))
        axi.plot([x, x], [2.38, 2.82], color="#1f5f8b", lw=1)
    axi.text(5.0, 5.6, "Cu", ha="center", fontsize=9, color="k")
    axi.text(5.0, 3.0, "乙醇 ＋ 双极化子", ha="center", fontsize=8,
             color="#1f5f8b", va="bottom")
    axi.annotate("", xy=(5, 5.2), xytext=(5, 4.5),
                 arrowprops=dict(arrowstyle="-|>", color="crimson", lw=1))
    axi.set_title("界面放大：垂直电流穿过乙醇", fontsize=8)
    fig.savefig(os.path.join(ROOT, "data", "fig_stack_2d.png"), dpi=160)
    plt.close(fig)


def fig_3d():
    fig = plt.figure(figsize=(10, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_axis_off()

    def plate(x, size=3.0, cu=0.06, ins=0.38):
        # 竖立：厚度沿 x，面在 yz 平面
        ax.bar3d(x, 0, 0, ins, size, size, color=INS, edgecolor="k",
                 linewidth=0.5, alpha=0.9)
        ax.bar3d(x - cu, 0, 0, cu, size, size, color=CU, edgecolor="k",
                 linewidth=0.3, alpha=0.95)
        ax.bar3d(x + ins, 0, 0, cu, size, size, color=CU, edgecolor="k",
                 linewidth=0.3, alpha=0.95)
        # 顶边缘连通
        ax.bar3d(x - cu, 0, size, ins + 2 * cu, size, 0.10, color=CU,
                 edgecolor="k", linewidth=0.3, alpha=0.95)

    xs = [0.0, 1.6, 3.2]
    for x in xs:
        plate(x)

    # 垫片
    for y, z in [(0.3, 0.3), (2.7, 0.3), (0.3, 2.7), (2.7, 2.7)]:
        for x in xs[:-1]:
            ax.scatter([x + 0.87], [y], [z], color="#7a7a7a", s=22)

    # 乙醇（透明，直接润湿两侧铜面，无空隙）
    for x in xs[:-1]:
        gx = x + 0.44          # 左铜面外沿
        gw = 1.10              # 填满到右铜面
        ax.bar3d(gx, 0, 0, gw, 3.0, 3.0, color=MEOH, alpha=0.25,
                 edgecolor="none")
        # 绝缘垫片框（四条边，加厚，兼密封）
        b = 0.25
        for (gy, gz, dy, dz) in [(0, 0, 3.0, b), (0, 3.0 - b, 3.0, b),
                                 (0, 0, b, 3.0), (3.0 - b, 0, b, 3.0)]:
            ax.bar3d(gx, gy, gz, gw, dy, dz,
                     color="#9a9a9a", alpha=0.55, edgecolor="none")

    # 电流（水平，画在顶部便于看见）
    ax.quiver(-0.8, 1.5, 3.55, 4.9, 0, 0, color="crimson",
              arrow_length_ratio=0.06, lw=1.8)
    ax.text(-1.3, 1.5, 3.75, "I", color="crimson", fontsize=12, fontweight="bold")

    ax.text(0.2, 3.6, 2.0, "乙醇间隙", color="#1f5f8b", fontsize=10)
    ax.text(1.2, -0.7, 0.0, "Cu 镀层", color=CU, fontsize=10)
    ax.text(-0.8, 0.0, 3.85, "顶边连通", color=CU, fontsize=9)
    ax.text(3.6, 1.0, 3.9, "叠压 →", color="crimson", fontsize=10)
    ax.set_title("镀铜铝板叠层：三维（竖立，N 片、N−1 个乙醇间隙）", fontsize=12)
    ax.view_init(elev=18, azim=-62)
    fig.savefig(os.path.join(ROOT, "data", "fig_stack_3d.png"), dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig_2d()
    fig_3d()
    print("saved data/fig_stack_2d.png, data/fig_stack_3d.png")

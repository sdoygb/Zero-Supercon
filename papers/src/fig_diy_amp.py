#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_diy_amp.py —— DIY 跨阻放大器（GΩ 测量）＋ ADC/单片机 框图"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams["font.sans-serif"] = ["Arial Unicode MS"]
plt.rcParams["axes.unicode_minus"] = False


def main():
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # 样品 R_x
    ax.add_patch(Rectangle((1.0, 4.6), 2.2, 1.2, fc="#e8e8e8", ec="k", lw=1.2))
    ax.text(2.1, 5.2, "样品 $R_x$\n(叠层)", ha="center", va="center", fontsize=9)
    # V_s 源
    ax.add_patch(Circle((1.0, 3.2), 0.45, fc="white", ec="k", lw=1.2))
    ax.text(1.0, 3.2, "$V_s$", ha="center", va="center", fontsize=10)
    ax.plot([1.0, 1.0], [3.65, 4.6], "k", lw=1.2)
    ax.plot([1.0, 3.6], [3.2, 3.2], "k", lw=1.2)
    ax.plot([3.6, 3.6], [3.2, 4.0], "k", lw=1.2)

    # 运放（三角）
    op = Polygon([[4.2, 3.2], [4.2, 5.2], [5.6, 4.2]], closed=True,
                 fc="#f7f7f7", ec="k", lw=1.4)
    ax.add_patch(op)
    ax.text(4.55, 4.3, "−", fontsize=12, ha="center")
    ax.text(4.55, 4.75, "+", fontsize=12, ha="center")
    # 输入
    ax.plot([3.6, 4.2], [4.0, 4.0], "k", lw=1.2)   # 反相输入
    ax.plot([4.2, 4.2], [4.75, 4.75], "k", lw=0)   # 占位
    ax.plot([3.6, 3.6], [4.0, 4.75], "k", lw=1.2)
    ax.plot([3.6, 4.2], [4.75, 4.75], "k", lw=1.2)  # 同相输入接地
    ax.plot([2.0, 2.0], [4.75, 5.6], "k", lw=1.0, ls=":")
    ax.text(1.6, 5.55, "GND", fontsize=8)

    # 反馈 R_f / C_f
    ax.plot([4.2, 4.2], [4.0, 6.2], "k", lw=1.2)
    ax.plot([5.6, 5.6], [4.2, 6.2], "k", lw=1.2)
    ax.plot([4.2, 4.9], [6.2, 6.2], "k", lw=1.2)
    ax.plot([5.0, 5.6], [6.2, 6.2], "k", lw=1.2)
    ax.add_patch(Rectangle((4.9, 6.0), 0.1, 0.4, fc="white", ec="k", lw=1.0))
    ax.text(5.35, 6.45, "$R_f=1$ GΩ（可切换）", ha="center", fontsize=9)
    ax.text(5.35, 5.75, "$C_f=1$ pF", ha="center", fontsize=8)

    # 输出 → ADC
    ax.plot([5.6, 6.6], [4.2, 4.2], "k", lw=1.2)
    ax.add_patch(Rectangle((6.6, 3.6), 1.5, 1.2, fc="#c6dbef", ec="k", lw=1.2))
    ax.text(7.35, 4.2, "ADC\nADS1115\n16 bit", ha="center", va="center", fontsize=8)
    ax.plot([8.1, 8.6], [4.2, 4.2], "k", lw=1.2)
    ax.add_patch(Rectangle((8.6, 3.4), 1.6, 1.6, fc="#a1d99b", ec="k", lw=1.2))
    ax.text(9.4, 4.2, "MCU\nESP32\n(I²C)", ha="center", va="center", fontsize=8)

    # 显示/SD
    ax.add_patch(Rectangle((10.5, 4.4), 1.0, 0.7, fc="#fee0b6", ec="k", lw=1.0))
    ax.text(11.0, 4.75, "OLED", ha="center", va="center", fontsize=7)
    ax.add_patch(Rectangle((10.5, 3.3), 1.0, 0.7, fc="#fee0b6", ec="k", lw=1.0))
    ax.text(11.0, 3.65, "SD 卡", ha="center", va="center", fontsize=7)
    ax.plot([10.2, 10.5], [4.75, 4.75], "k", lw=1.0)
    ax.plot([10.2, 10.5], [3.65, 3.65], "k", lw=1.0)

    # 温控
    ax.add_patch(Rectangle((8.6, 1.6), 1.6, 1.2, fc="#f4a582", ec="k", lw=1.2))
    ax.text(9.4, 2.2, "Peltier\nPID/PWM", ha="center", va="center", fontsize=8)
    ax.plot([9.4, 9.4], [3.4, 2.8], "k", lw=1.0)

    # 公式
    ax.text(6.0, 1.0, r"$V_o = V_s \cdot R_f / R_x$", ha="center", fontsize=13)
    ax.text(6.0, 0.45, "例：$R_x=1.45$ GΩ、$R_f=1$ GΩ、$V_s=2.5$ V → $V_o\\approx1.72$ V（ADS1115 直接读）",
            ha="center", fontsize=9)

    # 保护环提示
    ax.text(2.1, 7.3, "保护环（guard）绕输入走线；聚四氟乙烯支架；清洗助焊剂",
            ha="center", fontsize=8, color="#555")
    ax.set_title("DIY GΩ 测量：跨阻放大器 ＋ 16 位 ADC ＋ 单片机", fontsize=13)
    fig.savefig(os.path.join(ROOT, "data", "fig_diy_amp.png"), dpi=160)
    plt.close(fig)
    print("saved data/fig_diy_amp.png")


if __name__ == "__main__":
    main()

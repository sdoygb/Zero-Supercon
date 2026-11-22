#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_geometry.py —— 界面超导的几何示意：面内电流、跨面配对

产出 data/fig_geometry.png：
  左：A/B/A/B 叠层，电流沿面内（x），配对跨界面（z），界面模偏振沿 z；
  右：两种四端测量——面内（测零阻）与跨面（测 Josephson IcRn）。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

HF = [0.6, 1.0, 0.6, 1.0]              # 从下到上：A B A B
COL = ["#9ecae1", "#f4a582", "#9ecae1", "#f4a582"]
NAMES = ["A（轻网络）", "B（重金属）", "A（轻网络）", "B（重金属）"]


def draw_stack(ax, x0, w, labels=True):
    y = 0.0
    ys = []
    for h, c, nm in zip(HF, COL, NAMES):
        ax.add_patch(Rectangle((x0, y), w, h, facecolor=c, edgecolor="k", lw=1.2, alpha=0.9))
        if labels:
            ax.text(x0 - 0.08, y + h / 2, nm, ha="right", va="center", fontsize=9)
        ys.append((y, y + h))
        y += h
    for i in range(len(HF) - 1):
        yy = ys[i][1]
        ax.plot([x0, x0 + w], [yy, yy], color="k", lw=2.2)
    return y, ys


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=160)

# 左：几何
ax = axes[0]
x0, w = 0.0, 3.4
ytop, ys = draw_stack(ax, x0, w)
ax.add_patch(FancyArrowPatch((x0 + 0.15, ytop + 0.42), (x0 + w - 0.15, ytop + 0.42),
                             arrowstyle="-|>", mutation_scale=20, lw=2.2, color="#c0392b"))
ax.text(x0 + w / 2, ytop + 0.62, "电流：面内（平行界面）", ha="center", fontsize=9.5, color="#c0392b")
# 配对跨界面（A/B 之间的界面，y=1.6）
yi = ys[1][1]
xc = x0 + 2.45
ax.plot([xc, xc], [yi - 0.28, yi + 0.28], color="k", lw=1.2, ls=":")
ax.scatter([xc], [yi - 0.28], s=95, c="#2c3e50", zorder=5)
ax.scatter([xc], [yi + 0.28], s=95, c="#2c3e50", zorder=5)
ax.text(xc + 0.15, yi - 0.30, "配对：跨界面（z）", fontsize=9, va="center")
# 界面模偏振
xp = x0 + 0.55
ax.add_patch(FancyArrowPatch((xp, yi - 0.32), (xp, yi + 0.32),
                             arrowstyle="<|-|>", mutation_scale=14, lw=1.8, color="#16a085"))
ax.text(xp + 0.05, yi + 0.45, "界面模偏振（z）", fontsize=8.5, color="#16a085", ha="left")
ax.set_xlim(x0 - 1.55, x0 + w + 1.1)
ax.set_ylim(-0.25, ytop + 0.95)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("几何：电流沿界面，配对跨界面", fontsize=11)

# 右：两种测量
ax = axes[1]
ytop, ys = draw_stack(ax, x0, w, labels=False)
# 面内四端：沿上层的 B 内
yin = ys[3][0] + HF[3] / 2
ax.add_patch(FancyArrowPatch((x0 - 0.45, yin), (x0 + w + 0.05, yin),
                             arrowstyle="-|>", mutation_scale=16, lw=2.2, color="#c0392b"))
ax.text(x0 + w / 2, yin + 0.18, "面内四端：R=0、Ic", ha="center", fontsize=9, color="#c0392b")
# 跨面四端：竖箭头在右侧
ax.add_patch(FancyArrowPatch((x0 + w + 0.45, 0.15), (x0 + w + 0.45, ytop - 0.15),
                             arrowstyle="-|>", mutation_scale=16, lw=2.2, color="#2980b9"))
ax.text(x0 + w + 0.58, ytop / 2, "跨面四端：\nJosephson IcRn", fontsize=9, color="#2980b9", va="center")
ax.text(x0 + w / 2, -0.32, "（层数只是示意；真实界面须原子级接触）", ha="center", fontsize=8, color="0.35")
ax.set_xlim(x0 - 0.7, x0 + w + 1.9)
ax.set_ylim(-0.6, ytop + 0.5)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("两种测量：面内 vs 跨面", fontsize=11)

fig.tight_layout()
fig.savefig("data/fig_geometry.png", bbox_inches="tight")
print("saved data/fig_geometry.png")

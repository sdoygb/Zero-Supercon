#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig_window.py —— (lambda, omega_log) 平面上的常温窗口与 1058 条材料"""
import csv, json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"
OUT = os.path.join(ROOT, "data", "fig_rt_window.png")

plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

rows = list(csv.DictReader(open(SRC)))
mp = list(csv.DictReader(open("/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv")))
lam = np.array([float(r["lam"]) for r in rows])
w = np.array([float(r["wlog_K"]) for r in rows])
tc = np.array([float(r["Tc_K"]) if r["Tc_K"] else 0.0 for r in rows])
ok = (lam > 0.02) & (w > 1)

def wreq(l, ms=0.10, T=300.0):
    l = np.asarray(l, float)
    den = l - ms * (1 + 0.62 * l)
    out = np.full_like(l, np.nan)
    m = den > 1e-6
    out[m] = 1.2 * T * np.exp(1.04 * (1 + l[m]) / den[m])
    return out

ml = np.array([float(r["lam"]) for r in mp if r["lam"] and r["wlog_K"]])
mw = np.array([float(r["wlog_K"]) for r in mp if r["lam"] and r["wlog_K"]])
mok = (ml > 0.02) & (mw > 1)

fig, ax = plt.subplots(figsize=(7.6, 5.6), dpi=160)
ax.scatter(ml[mok], mw[mok], s=3, c="0.55", alpha=0.25, linewidths=0, label="MP 7678（灰，λ 为模型值）")
sc = ax.scatter(lam[ok], w[ok], c=np.log10(np.maximum(tc[ok], 0.05)), s=9,
                cmap="viridis", alpha=0.85, linewidths=0)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.02, 30); ax.set_ylim(4, 2e4)
ax.set_xlabel(r"电声耦合 $\lambda$")
ax.set_ylabel(r"$\omega_{\log}$ (K)")

xs = np.logspace(-1.5, 1.5, 400)
ys = wreq(xs, 0.10)
ax.plot(xs, ys, "r-", lw=1.8, label=r"$T_c=300$ K（$\mu^*=0.10$）")
ax.fill_between(xs, ys, 2e4, color="red", alpha=0.07)
ax.axhline(300 / 0.275, color="r", ls=":", lw=1.1)
ax.text(0.025, 300 / 0.275 * 1.12, r"上限 $T_c<0.275\,\omega_{\log}$ $\Rightarrow$ $\omega_{\log}>1091$ K",
        color="r", fontsize=8.5)

for name, x, y, dx, dy in [("MoN", 1.367, 307.3, 5, -9),
                           (r"B$_2$Mg", 0.696, 846.3, 6, -11),
                           (r"B$_2$CN", 0.555, 912.4, -30, 4)]:
    ax.annotate(name, (x, y), textcoords="offset points", xytext=(dx, dy),
                fontsize=8, color="0.15")

cb = fig.colorbar(sc, ax=ax, pad=0.01)
cb.set_label(r"$\log_{10}(T_c/\mathrm{K})$")
ax.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
ax.set_title("常温窗口与常压材料（JARVIS 1058 ＋ MP 7678）", fontsize=11)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT)
print("saved", OUT)

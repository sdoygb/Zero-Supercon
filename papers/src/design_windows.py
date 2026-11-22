#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""design_windows.py —— 两种金属／交替叠层的能隙矩阵与设计窗口

输出 data/design_windows.json：
- two_band_delta：带边失配 δ 的带间核与 Tc（MgB2 型、金刚石型）
- mass_match：质量失配 → N_eff → Tc
- mu_star：多子带库仑赝势矩阵 λ_eff
- scenarios：端到端情景表
- bipolaron_3d / bipolaron_2d：双极化子 BEC／BKT 温度
- cu_overlayer：Cu ＋ 高 ω 覆盖层的两带 Tc

所有温度 K，能量 eV。
"""
from __future__ import annotations
import json
import math
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "design_windows.json")

HBAR = 1.054571817e-34
KB = 1.380649e-23
EV = 1.602176634e-19
ME = 9.1093837e-31
MS0 = 0.12      # 带内库仑赝势
MS12 = 0.04     # 带间库仑赝势


def ad(lam, w=1916.0, ms=MS0):
    """Allen–Dynes 形式的 Tc。"""
    den = lam - ms * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    return w / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def lam_req(w, T=300.0, ms=MS0):
    """给定 ω 与目标 Tc，反解 λ；天花板不足时返回 None。"""
    if 0.275 * w < T:
        return None
    rhs = math.log(w * (1 / 1.2) / T)
    return (1.04 + rhs * ms) / (rhs * (1 - 0.62 * ms) - 1.04)


def L(T, w):
    return math.log(1.13 * w / T) if 0 < T < 1.13 * w else -1e9


def L12(T, w, d):
    """带边失配 δ 下的带间配对核。"""
    if d <= 0:
        return L(T, w)
    f = lambda x: math.tanh(x / (2 * T)) / (x + d / 2)
    return quad(f, 0, w, limit=200)[0]


def two_band_tc(l11, w1, l22, w2, l12, w12, ms0=MS0, ms12=MS12, delta=0.0):
    """两带能隙矩阵的最大本征值跨过 1 的温度（取最高根）。"""
    def eig(T):
        M = np.array([
            [(l11 - ms0) * L(T, w1), (l12 - ms12) * L12(T, w12, delta)],
            [(l12 - ms12) * L12(T, w12, delta), (l22 - ms0) * L(T, w2)],
        ])
        return np.linalg.eigvalsh(M)[-1]
    Ts = np.logspace(math.log10(0.2), math.log10(max(w1, w2) * 1.1), 800)
    roots = [Ts[i] for i in range(1, len(Ts))
             if (eig(Ts[i]) - 1) * (eig(Ts[i - 1]) - 1) < 0]
    return max(roots) if roots else 0.0


def neff(eta, wD=0.165, D=0.0056, cap=40.0):
    """质量失配 η 给的有效子带数。"""
    if eta <= 0:
        return cap
    return min(cap, math.sqrt(wD / (D * eta)))


def t_bec3d(n_cm3, mstar):
    n = n_cm3 * 1e6
    return 3.31 * HBAR ** 2 * n ** (2 / 3) / (mstar * ME * KB)


def t_bkt2d(n2d_cm2, mstar):
    n = n2d_cm2 * 1e4
    TF = 2 * math.pi * HBAR ** 2 * n / (mstar * ME * KB)
    T = TF / 10
    for _ in range(200):
        T = TF / math.log(380 * TF / T)
    return T


def main():
    out = {}

    # 1) 带边失配
    for name, w, l0, l12 in (("mgb2", 846.0, 0.30, 0.20),
                             ("diamond", 1916.0, 0.30, 0.20)):
        rows = []
        for frac in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8):
            d = frac * w
            rows.append({"delta_over_w": frac, "delta_K": d,
                         "Tc": two_band_tc(l0, w, 0.1, w, l12, w, delta=d)})
        out["two_band_delta_" + name] = rows

    # 2) 质量失配
    mm = []
    for eta in (0.0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.20):
        n = neff(eta)
        lam = (0.46 - MS0) + (n - 1) * (0.05 - MS12) * 0.6
        mm.append({"eta": eta, "Neff": n, "Tc": ad(lam)})
    out["mass_match"] = mm

    # 3) μ* 矩阵
    mu = []
    for l12 in (0.05, 0.10, 0.15, 0.20):
        row = {"lam12": l12}
        for m12 in (0.0, 0.02, 0.04, 0.06, 0.08):
            lam = (0.46 - MS0) + 39 * (l12 - m12)
            row["mu%.2f" % m12] = ad(lam)
        mu.append(row)
    out["mu_star"] = mu

    # 4) 端到端情景
    scen = []
    for name, l12, m12, eta, p in (
            ("best", 0.10, 0.04, 0.01, 1.0),
            ("optimistic", 0.08, 0.05, 0.01, 0.8),
            ("nominal", 0.06, 0.05, 0.03, 0.6),
            ("pessimistic", 0.05, 0.06, 0.10, 0.5),
            ("worst", 0.04, 0.07, 0.10, 0.5)):
        n = neff(eta)
        lam = (0.46 - MS0) + (n - 1) * (l12 - m12) * p
        scen.append({"name": name, "lam12": l12, "mu12": m12,
                     "eta": eta, "Neff": n, "parity": p,
                     "lam_tot": lam, "Tc": ad(lam)})
    out["scenarios"] = scen

    # 5) 双极化子 BEC / BKT
    out["bipolaron_3d"] = [
        {"n_cm3": n, **{f"m{m}": t_bec3d(n, m) for m in (2, 5, 10, 20)}}
        for n in (1e20, 3e20, 1e21, 3e21, 1e22)]
    out["bipolaron_2d"] = [
        {"n2d_cm2": n, **{f"m{m}": t_bkt2d(n, m) for m in (2, 3, 5, 10)}}
        for n in (1e13, 3e13, 1e14, 3e14, 1e15)]

    # 6) Cu + 高 ω 覆盖层
    cu = []
    for l12 in (0.1, 0.2, 0.3, 0.5):
        cu.append({"lam12": l12,
                   "Tc_diamond": two_band_tc(0.158, 223.0, 0.1, 1916.0, l12, 1916.0),
                   "Tc_aln": two_band_tc(0.158, 223.0, 0.1, 900.0, l12, 900.0)})
    out["cu_overlayer"] = cu

    # 7) 关键门槛
    out["thresholds"] = {
        "lam_req_1916": lam_req(1916.0),
        "lam_req_5180": lam_req(5180.0),
        "lam_req_846": lam_req(846.0),
        "lam_req_223": lam_req(223.0),
    }

    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)
    for k in ("thresholds",):
        print(k, out[k])
    print("best scenario:", [s for s in out["scenarios"] if s["name"] == "best"][0])


if __name__ == "__main__":
    main()

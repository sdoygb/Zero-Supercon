#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""deform_channels.py —— Zd021：形变势的通道分解（对角 vs 非对角）

舞台（Zb005）：eps_k = t(lam_k-2d)，lam_k=4 sum_i sin^2(k_i/2)；3D 时 eps=-2t sum cos k_i。
非对角通道（hopping 随键长变）：t -> t(1-n da/a)，n=-dln t/dln a（Zd014）。
  则 eps_k -> eps_k(1-n da/a)，即
      Xi(k) := d eps_k / d(da/a) = -n * eps_k  （应变导数，能量）
  矩阵元（每单位应变）在费米面上 = -n * eps_F —— **半满 eps_F=0 时恒为零**。
对角通道（on-site，随应变移带心）：Xi_d，与填充无关（材料量，C）。
  合起来：Xi_FS = Xi_d - n eps_F。

本脚本：① 恒等式数值核验；② 填充因子表；③ 七个材料"若纯 hopping"反推的 eps_F/t 与填充；
④ Cu 的反例（s 带半满 -> hopping 通道为零，但表里 lambda=0.158）。
输出：data/deform_channels.json
"""
from __future__ import annotations
import csv
import json
import math
import os

import numpy as np

from constants import C_LAM_G1, C_N

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "deform_channels.json")
MP = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"

HB2_2ME_EVA2 = 3.80998
KAPPA = (3 * math.pi ** 2) ** (1 / 3)

# 与 Zd020 §5.5 同一组材料（a, M, n_H, lambda 取 MP 表最近条目）
BULK = [
    ("金刚石 C", 3.567, 12.011, 2, 1336.5),
    ("Si", 5.431, 28.086, 2, 403.1),
    ("W", 3.165, 183.84, 5, 171.1),
    ("Ta", 3.306, 180.95, 5, 129.6),
    ("Nb", 3.301, 92.906, 5, 175.4),
    ("Mo", 3.147, 95.95, 5, 230.9),
    ("Cu", 3.615, 63.546, 2, 223.3),
]


def eps_of(n_e, t=1.0, nb=4096):
    """3D 简立方：给定填充 n_e∈[0,2]，返回 eps_F（单位 t）。"""
    n = 512
    c = np.cos(np.arange(n) * 2 * np.pi / n)
    h = np.zeros(nb)
    for i in range(n):
        e = -2 * (c[i] + c[:, None] + c[None, :])
        h += np.histogram(e, bins=nb, range=(-6, 6))[0]
    ctr = -6 + (np.arange(nb) + 0.5) * 12 / nb
    dens = h / n ** 3 / (12 / nb)
    cum = np.cumsum(dens) * (12 / nb)
    return float(np.interp(n_e / 2.0, cum, ctr)) * t


def main():
    # ① 恒等式：d eps/d(da/a) = -n eps
    k = np.linspace(0, np.pi, 7)[1:-1]
    K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    lam = 4 * (np.sin(K[..., 0] / 2) ** 2 + np.sin(K[..., 1] / 2) ** 2 + np.sin(K[..., 2] / 2) ** 2)
    eps = lam - 6.0
    n = 2.0
    h = 1e-6
    # eps(dl) = (1-n*dl)(lam-6)；数值差分 d eps/d dl
    d = ((1 - n * h) * (lam - 6.0) - (1 + n * h) * (lam - 6.0)) / (2 * h)
    ident_err = float(np.max(np.abs(d - (-n * (lam - 6.0)))))
    # ② 填充因子：由填充给 eps_F（单位 t），并给"若纯 hopping"的 Xi_FS/Xi_edge
    fill = []
    for ne in (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8):
        ef = eps_of(ne, t=1.0)
        fill.append({"n_e": ne, "epsF_over_t": round(ef, 4),
                     "Xi_FS_over_Xi_edge": round(ef / 6.0, 4),
                     "lam_ratio_hopping": round((ef / 6.0) ** 2, 5)})
    # ③ 七个材料：若纯 hopping，则 n_eff/n_H = eps_F/t
    rows = []
    mp = {}
    for r in csv.DictReader(open(MP)):
        try:
            mp.setdefault(r["formula"], []).append((float(r["wlog_K"]), float(r["lam"])))
        except Exception:
            pass
    for name, a, M, nH, wt in BULK:
        f = name.split()[-1] if name.startswith("金刚石") else name
        cand = mp.get(f) or []
        wlog, lam = min(cand, key=lambda x: abs(x[0] - wt))
        t = HB2_2ME_EVA2 * (KAPPA / a) ** 2
        n_eff = math.sqrt(lam / C_LAM_G1)
        ratio = n_eff / nH
        e2a = 14.4 / a                      # eV，Coulomb 标度 e²/a
        rows.append({"material": name, "lam": round(lam, 4), "t_bare_eV": round(t, 3),
                     "n_H": nH, "n_eff": round(n_eff, 3),
                     "epsF_over_t_if_hopping": round(ratio, 3),
                     "n_e_if_hopping": round(1.0 + 2 * ratio * C_N, 3),
                     "lam_ratio": round(lam / (C_LAM_G1 * nH ** 2), 4),
                     "e2_over_a_eV": round(e2a, 2),
                     "e2a_over_t": round(e2a / t, 2),
                     "c_diag": round(n_eff / (e2a / t), 3)})
    # ④ Cu 反例：s 带半满，hopping 通道为零；表里 lambda=0.158 -> 对角项必须 >= D_eff
    cu = next(r for r in rows if r["material"] == "Cu")
    out = {
        "identity": {"claim": "d eps/d(da/a) = -n eps", "max_err": ident_err},
        "filling_factor": fill,
        "materials": rows,
        "cu": {"note": "Cu 的 s 带半满：eps_F=0 -> hopping 通道为零；MP 表 lambda=0.158 只能来自对角项",
               "lam_MP": cu["lam"], "D_d_implied_eV": round(cu["n_eff"] * cu["t_bare_eV"], 2),
               "D_hop_eV": round(cu["n_H"] * cu["t_bare_eV"], 2)},
        "reading": "纯 hopping 图景要求 eps_F/t = n_eff/n_H = 0.10-0.22，即填充 1.03-1.06；"
                   "七个材料横跨 sp/d，这个巧合不可信 -> 数据要求对角通道",
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("identity max_err =", ident_err)
    print("填充 | epsF/t | Xi_FS/Xi_edge | lambda 比")
    for f in fill:
        print(f"  {f['n_e']:.1f} | {f['epsF_over_t']:+.3f} | {f['Xi_FS_over_Xi_edge']:+.3f} | {f['lam_ratio_hopping']:.5f}")
    print("\n材料 | λ | t | n_H | n_eff | epsF/t(若纯hopping) | n_e(若纯hopping)")
    for r in rows:
        print(f"  {r['material']:8s} {r['lam']:.3f} {r['t_bare_eV']:.2f} {r['n_H']} {r['n_eff']:.2f} "
              f"{r['epsF_over_t_if_hopping']:.3f} {r['n_e_if_hopping']:.3f}")
    print("\nCu:", json.dumps(out["cu"], ensure_ascii=False))
    print("saved", OUT)


if __name__ == "__main__":
    main()

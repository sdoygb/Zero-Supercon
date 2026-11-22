#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""two_component_parity.py —— 双组份（双带）宇称平衡的候选扫描

模型：两个带
  a（轻网络）：ω_a（网络的光学模）、N_a（界面态，C）
  b（金属）：  ω_b（金属声子）、N_b（Zd026 的 d 带 DOS）
  λ_aa、λ_bb（带内）、λ_ab=λ_ba（带间，经界面键 D_I）
宇称平衡：带间耦合的有效部分 λ_ab^eff = λ_ab·(2 w_A1)（Zd027 的判据 w_A1>1/2）
Tc：线性化双带方程的 2×2 行列式
  (λ_aa L_a −1)(λ_bb L_b −1) = λ_ab² L_a L_b，L_i=ln(1.13 ω_i/Tc)

输出：data/two_component_parity.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "two_component_parity.json")

K = 0.1660
MU = 0.12
# 金属（Zd026 / Zd029 的 N_B）
METALS = {"Hf": 0.67, "Ta": 0.92, "W": 1.54, "Nb": 2.05, "Mo": 4.66,
          "Ni": 1.63, "Co": 1.47, "Pt": 1.31, "Cu": 0.107, "Ag": 0.107}
# 轻网络：光学模（K）与界面宇称权重 w_A1（Zd027：sp³ 1/4、sp² 1/3、s 主导 >1/2）
NETWORKS = {
    "金刚石 C–C（sp³，1916 K）": (1916.0, 0.25),
    "sp² C–N（C₃N₄，1500 K）": (1500.0, 1 / 3),
    "hBN（B–N，1966 K）": (1966.0, 1 / 3),
    "石墨烯 G（2274 K）": (2274.0, 1 / 3),
    "s 带金属侧（Cu/Ag，1000 K）": (1000.0, 0.60),
}
# C–M 的 D_I（Zd029）
DI = {"Hf": 12.14, "Ta": 9.83, "W": 8.86, "Nb": 9.73, "Mo": 9.24,
      "Ni": 7.64, "Co": 8.25, "Pt": 8.83, "Cu": 7.53, "Ag": 8.67}
W_B = 350.0        # 金属声子（K），C
N_A = 0.30         # 轻网络侧的界面态 DOS（C）


def tc_two_band(w_a, l_aa, w_b, l_bb, l_ab, mu):
    """解最大本征值 =1 的 T：M=[[λaa La, λab La],[λab Lb, λbb Lb]]，La=ln(1.13ωa/T)。"""
    lab = max(l_ab - mu, 0.0)

    def eig(T):
        La = math.log(1.13 * w_a / T) if T < 1.13 * w_a else 0.0
        Lb = math.log(1.13 * w_b / T) if T < 1.13 * w_b else 0.0
        M = [[l_aa * La, lab * La], [lab * Lb, l_bb * Lb]]
        tr = M[0][0] + M[1][1]
        dt = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        return 0.5 * (tr + math.sqrt(max(tr * tr - 4 * dt, 0.0)))

    lo, hi = 0.5, 1.13 * w_a
    if eig(lo) < 1.0:
        return 0.0
    if eig(hi) > 1.0:
        return hi
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if eig(mid) > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    rows = []
    for net, (w_a, wA1) in NETWORKS.items():
        for m, nb in METALS.items():
            lam_ab = K * math.sqrt(N_A * nb) * DI[m] ** 2
            lam_ab_eff = lam_ab * (2 * wA1)          # 宇称平衡：只留偶分量
            lam_bb = K * nb * DI[m] ** 2
            lam_aa = K * N_A * DI[m] ** 2 * 0.1      # 轻网络侧带内（弱）
            tc = tc_two_band(w_a, lam_aa, W_B, lam_bb, lam_ab_eff, MU)
            rows.append({"network": net, "metal": m, "w_a_K": w_a, "w_A1": round(wA1, 2),
                         "N_b": nb, "D_I": DI[m], "lam_ab": round(lam_ab, 2),
                         "lam_ab_eff": round(lam_ab_eff, 2), "lam_bb": round(lam_bb, 2),
                         "Tc": round(tc), "parity_ok": wA1 > 0.5})
    rows.sort(key=lambda r: -r["Tc"])
    print("Top 12（按 Tc）：")
    print("  %-22s %-4s %6s %5s %7s %7s %6s %s" % ("网络", "金属", "ω_a", "w_A1", "λab_eff", "λbb", "Tc", "宇称?"))
    for r in rows[:12]:
        print("  %-22s %-4s %6.0f %5.2f %7.2f %7.2f %6d %s"
              % (r["network"], r["metal"], r["w_a_K"], r["w_A1"], r["lam_ab_eff"],
                 r["lam_bb"], r["Tc"], "过" if r["parity_ok"] else "不过"))
    json.dump({"rows": rows, "mu_star": MU, "K": K, "W_b": W_B, "N_a": N_A},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

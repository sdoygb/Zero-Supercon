#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""omega_valence.py —— 界面模频率的价加权（Zc038 命题 6.4）

命题 6.4：碳价闭合 sum_j A_ij + n_i = 4 ⟹ 一个碳参与 nu 个键。
故在"单键振子"里，碳的质量份额是 m_C/nu，不是整个原子。
Morse 曲率 K_I = 2 a² D_b；界面键 C–M 的约化质量 μ_I = (m_C/nu·m_M)/(m_C/nu + m_M)。

核验：用 Raman 实测（金刚石 1332 cm^-1、石墨烯 1580 cm^-1）与三种取法对照：
  全原子（旧）、价加权（本脚本）、实测。
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "omega_valence.json")

EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
HB = 6.582119569e-16          # eV·s
KB = 8.617333e-5              # eV/K
CM_K = 1.4387769              # 1 cm^-1 对应的 K


def omega_K(KI_Nm, mu_amu):
    w = math.sqrt(KI_Nm / (mu_amu * AMU))     # rad/s
    return HB * w / KB                         # K


def omega_cm(KI_Nm, mu_amu):
    return omega_K(KI_Nm, mu_amu) / CM_K


def main():
    mC, mAg = 12.0, 107.87
    nets = [
        # 名称, D_b(eV), a(Å^-1), nu, 实测 C–C 模(cm^-1)
        ("金刚石 sp3", 3.6, 2.0, 4, 1332.0),
        ("sp2 石墨烯", 4.2, 2.2, 3, 1580.0),
        ("C3N4", 3.0, 2.0, 3, None),
    ]
    rows = []
    print("%-11s %6s %6s %5s | %8s %8s %8s" % ("网络", "K_I", "ω全胞", "ω价加", "实测", "价加权比", "全原子比"))
    for name, Db, a, nu, exp in nets:
        KI = 2 * Db * a ** 2 * EV_A2
        m_full = mC
        m_part = mC / nu
        mu_full = m_full / 2                      # 键振子：两原子约化质量 = m/2
        mu_part = m_part / 2
        w_full = omega_cm(KI, mu_full)
        w_part = omega_cm(KI, mu_part)
        # 界面键（C–M，M=Ag）
        muI_full = mC * mAg / (mC + mAg)
        muI_part = m_part * mAg / (m_part + mAg)
        wI_full = omega_K(KI, muI_full)
        wI_part = omega_K(KI, muI_part)
        rows.append({"net": name, "D_b_eV": Db, "a_perA": a, "nu": nu, "K_I_Nm": round(KI),
                     "omega_full_atom_cm": round(w_full), "omega_valence_cm": round(w_part),
                     "omega_exp_cm": exp,
                     "omega_I_full_atom_K": round(wI_full), "omega_I_valence_K": round(wI_part),
                     "mu_I_valence_amu": round(muI_part, 2)})
        print("%-11s %6.0f %8.0f %8.0f %8s %8.2f %8.2f"
              % (name, KI, w_full, w_part, exp if exp else "—",
                 (w_part / exp if exp else float("nan")), (w_full / exp if exp else float("nan"))))
    print("\n界面键（C–Ag）的 ω_I：")
    for r in rows:
        print("  %-11s 全原子 %.0f K → 价加权 %.0f K（μ_I=%.2f amu）"
              % (r["net"], r["omega_I_full_atom_K"], r["omega_I_valence_K"], r["mu_I_valence_amu"]))
    json.dump({"CM_K": CM_K, "rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

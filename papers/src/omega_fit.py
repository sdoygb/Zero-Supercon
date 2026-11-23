#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""omega_fit.py —— 用 Raman 实测定"有效振子质量"

两种取法（全原子 / 价加权）都不是定理：
  全原子：μ = m_C/2（键振子，两原子各贡献整个原子）
  价加权：μ = (m_C/ν)/2（一个碳参与 ν 个键，每键质量份额 m_C/ν）
实测（C–C 光学模）：金刚石 1332 cm^-1、石墨烯 1580 cm^-1。

本脚本：①给两种取法的偏差；②反解"有效振子质量 μ_eff"（由实测与 K_I）；
③核验有效质量是否与 m_C/ν 或 m_C/2 同量级。
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "omega_fit.json")

EV_A2 = 1.602176634e-19 / 1e-20
AMU = 1.66054e-27
HB = 6.582119569e-16
KB = 8.617333e-5
CM_K = 1.4387769


def omega_cm(KI_Nm, mu_amu):
    w = math.sqrt(KI_Nm / (mu_amu * AMU))
    return (HB * w / KB) / CM_K


def mu_from_omega(KI_Nm, omega_cm_obs):
    """反解 μ（amu）：ω_rad = 2πc·ω_cm"""
    w_dev = (omega_cm_obs * CM_K) * KB / HB            # rad/s
    return KI_Nm / (w_dev ** 2) / AMU


def main():
    mC, mAg = 12.0, 107.87
    nets = [("金刚石 sp3", 3.6, 2.0, 4, 1332.0), ("sp2 石墨烯", 4.2, 2.2, 3, 1580.0)]
    rows = []
    print("%-11s %5s %7s %8s %8s %8s %8s" % ("网络", "ν", "K_I", "ω实测", "μ_eff", "μ_C/2", "μ_(C/ν)/2"))
    for name, Db, a, nu, exp in nets:
        KI = 2 * Db * a ** 2 * EV_A2
        w_full = omega_cm(KI, mC / 2)
        w_part = omega_cm(KI, (mC / nu) / 2)
        mu_eff = mu_from_omega(KI, exp)
        rows.append({"net": name, "nu": nu, "K_I_Nm": round(KI), "omega_exp_cm": exp,
                     "omega_full_atom_cm": round(w_full), "omega_valence_cm": round(w_part),
                     "mu_eff_amu": round(mu_eff, 2),
                     "mu_full_amu": mC / 2, "mu_valence_amu": round((mC / nu) / 2, 2),
                     "dev_full": round(w_full / exp, 2), "dev_valence": round(w_part / exp, 2)})
        print("%-11s %5d %7.0f %8.0f %8.2f %8.1f %8.2f"
              % (name, nu, KI, exp, mu_eff, mC / 2, (mC / nu) / 2))
    print("\n判读：μ_eff 落在 (μ_full, μ_valence) 之间 ⟹ 两取法定的是上下界。")
    for r in rows:
        lo, hi = min(r["mu_full_amu"], r["mu_valence_amu"]), max(r["mu_full_amu"], r["mu_valence_amu"])
        inside = lo <= r["mu_eff_amu"] <= hi
        print("  %-11s μ_eff=%.2f ∈ [%.1f, %.1f] ? %s" % (r["net"], r["mu_eff_amu"], lo, hi, inside))
    json.dump({"rows": rows}, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("saved", OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lifetime.py —— Zd006 的寿命条件

tau = tau0 * exp(Ea/kT),  tau0 = 2*pi*hbar/omega_max
Ea_min = kT * ln(tau_op/tau0)
输出：data/lifetime.json
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "lifetime.json")
HBAR = 6.582119569e-16  # eV s
KB = 8.617333e-5        # eV/K
T = 300.0
KBT = KB * T


def tau0_of(wmax_ev):
    return 2 * math.pi * HBAR / wmax_ev


def ea_min(tau_op, wmax_ev):
    return KBT * math.log(tau_op / tau0_of(wmax_ev))


def main():
    ops = {"1 s": 1.0, "1 day": 86400.0, "1 yr": 3.156e7, "100 yr": 3.156e9}
    rows = []
    for name, tau in ops.items():
        rows.append({"tau_op": name, "Ea_min_eV": round(ea_min(tau, 0.2), 3)})
    wrows = []
    for w in (0.1, 0.2, 0.4):
        wrows.append({"omega_max_eV": w, "tau0_s": tau0_of(w),
                      "Ea_min_1yr_eV": round(ea_min(3.156e7, w), 3)})
    out = {"kT_eV": KBT, "rows": rows, "omega_scan": wrows}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("kT =", round(KBT, 5), "eV")
    for r in rows:
        print(f"  {r['tau_op']:7s} Ea_min = {r['Ea_min_eV']:.3f} eV")
    print("omega scan (1 yr):")
    for r in wrows:
        print(f"  omega={r['omega_max_eV']:.1f} eV  tau0={r['tau0_s']:.3e} s  Ea_min={r['Ea_min_1yr_eV']:.3f} eV")
    print("saved", OUT)


if __name__ == "__main__":
    main()

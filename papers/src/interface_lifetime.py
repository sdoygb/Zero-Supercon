#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_lifetime.py —— Zd017：网络保留界面的寿命

tau = tau0 exp(Ea/kT), tau0 = 2 pi hbar / w_max（Zd006）。
Ea 取 0.5–4.0 eV 的扫描；同时给"1 年"的阈值。
输出：data/interface_lifetime.json
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_lifetime.json")
HBAR = 6.582119569e-16
KB = 8.617333e-5
T = 300.0
WMAX = 0.2  # eV
TAU0 = 2 * math.pi * HBAR / WMAX


def tau(ea):
    return TAU0 * math.exp(ea / (KB * T))


def main():
    rows = [{"Ea_eV": ea, "tau_s": tau(ea),
             "log10_tau": round(math.log10(tau(ea)), 2)} for ea in (0.5, 0.8, 1.0, 1.26, 1.5, 2.0, 2.5, 3.0, 3.6)]
    one_yr = KB * T * math.log(3.156e7 / TAU0)
    out = {"tau0_s": TAU0, "kT_eV": KB * T, "Ea_1yr_eV": round(one_yr, 3), "rows": rows}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"tau0 = {TAU0:.3e} s; kT = {KB*T:.5f} eV; Ea(1 yr) = {one_yr:.3f} eV")
    for r in rows:
        print(f"  Ea={r['Ea_eV']:.2f} eV  tau={r['tau_s']:.2e} s  log10={r['log10_tau']}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

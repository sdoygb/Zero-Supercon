#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""superlattice_geometry.py —— Zd028：超晶格的三角约束（层厚 d、相干长度 ξ、界面密度）

三个量：
  ξ = hbar v_F/(π Δ)，Δ = 1.76 k_B T_c          —— 凝聚体的厚度
  f = min(1, 2ξ/d)                              —— 界面增强覆盖的体积分数
  n_int = 2/d                                   —— 单位长度的界面数（每层两个界面）
判据：体超导（非 S-I-S）要 f≈1；同时化学要 d ≥ d_chem（层不被界面吃掉）。
窗口存在 ⟺ d_chem < 2ξ。

输出：data/superlattice_geometry.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "superlattice_geometry.json")

HBAR = 1.0546e-34
E_CHARGE = 1.602176634e-19
KB_EV = 8.617333e-5
A_LAYER = 0.30e-9      # 一个原子层的厚度（Å 量级）


def xi_nm(Tc, vF):
    D = 1.76 * KB_EV * Tc                      # eV
    return HBAR * vF / (math.pi * D * E_CHARGE) * 1e9


def main():
    vFs = {"v_F=4.6e5": 4.6e5, "v_F=1.33e6": 1.33e6}
    tcs = (10, 50, 100, 300)
    xi = {k: {str(T): round(xi_nm(T, v), 1) for T in tcs} for k, v in vFs.items()}

    # 体积分数 f(d) 与界面密度，取 ξ=6.1 nm
    xi0 = xi_nm(300, 1.33e6)
    ds = (1, 2, 3, 5, 8, 12, 20, 50)
    frac = [{"d_nm": d, "layers": round(d * 1e-9 / A_LAYER, 1),
             "f": round(min(1.0, 2 * xi0 / d), 2),
             "n_int_per_nm": round(2.0 / d, 3)} for d in ds]

    # 窗口：d ∈ [d_chem, 2ξ]
    window = {"xi_300K_nm": round(xi0, 1), "d_max_nm": round(2 * xi0, 1),
              "examples": {f"d_chem={d}": ("窗口非空" if d < 2 * xi0 else "窗口为空")
                           for d in (2, 5, 10, 12, 15)}}

    out = {"xi_nm": xi, "frac": frac, "window": window,
           "note": "f=2ξ/d 是上界（均匀增强）；界面粗糙度、双界面模杂化、微带隙未计"}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))
    print("saved", OUT)


if __name__ == "__main__":
    main()

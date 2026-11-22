#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""migdal_zone.py —— Zd016 §3.5：Migdal 条件的适用区

判据：Γ := ω_D/ε_F。顶点修正 ~Γ，故单圈 Eliashberg 要求 Γ≪1。
表里的 ω_D、ε_F 是材料输入（C）；本脚本只算比值并给分档。

输出：data/migdal_zone.json
"""
from __future__ import annotations
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "migdal_zone.json")

# (类, ω_D[eV], ε_F[eV], 备注)
TABLE = [
    ("金刚石（C）", 0.165, 10.0, "光支 1332 cm⁻¹；宽 sp 带"),
    ("H$_3$S", 0.15, 7.0, "高压氢化物"),
    ("MgB$_2$", 0.08, 1.5, "E$_{2g}$ 模"),
    ("平带材料", 0.10, 0.5, "窄带/平带"),
    ("分子晶体", 0.08, 0.3, "窄带"),
]


def main():
    rows = []
    for name, w, ef, note in TABLE:
        g = w / ef
        band = "单圈够" if g < 0.10 else ("顶点修正可辨" if g < 0.20 else "失效区")
        rows.append({"class": name, "omega_eV": w, "epsF_eV": ef,
                     "Gamma": round(g, 3), "band": band, "note": note})
    out = {"criterion": "Gamma = omega_D/eps_F; vertex ~ Gamma; one-loop needs Gamma<<1",
           "rows": rows,
           "reading": "轻网络 Gamma=0.017-0.053；窄带/平带 Gamma=0.20-0.27"}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    for r in rows:
        print(f"{r['class']:10s} w={r['omega_eV']:.3f} epsF={r['epsF_eV']:5.1f}  "
              f"Gamma={r['Gamma']:.3f}  {r['band']}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""parity_form_factor.py —— Zd027：界面配对的宇称判据与 C–M 的形状因子

设置：C(001)/M(001) 界面，点群 C4v，宇称算符是镜面 M_z（z -> -z）。
轨道宇称（M_z）：
  s +, pz -, px/py +
  dz² +, dxz/dyz -, dxy/dx²-y² +
配对通道的宇称 = 两个轨道宇称之积；同轨道内（intra-band）恒为偶（偶×偶、奇×奇）。
界面 C 原子的悬挂键是 sp³ 杂化：(1/2)s + (√3/2)pz —— s 权重 1/4、pz 权重 3/4。

输出：data/parity_form_factor.json
"""
from __future__ import annotations
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "parity_form_factor.json")

PARITY = {"s": +1, "pz": -1, "px": +1, "py": +1,
          "dz2": +1, "dxz": -1, "dyz": -1, "dxy": +1, "dx2y2": +1}

# sp³ 杂化：沿 +z 的悬挂键
C_S = 0.5
C_PZ = math.sqrt(3) / 2


def channel(o1, o2):
    return PARITY[o1] * PARITY[o2]


def main():
    # 1) 轨道宇称表
    table = [{"orb": o, "parity_Mz": p} for o, p in PARITY.items()]

    # 2) 沿 z 的界面键：允许的轨道对（σ: m=0；π: m=±1）
    # 沿 z 的键只有这四类有重叠（其余由角因子为零）
    bonds = [("s", "dz2", "σ"), ("pz", "dz2", "σ"), ("px", "dxz", "π"), ("py", "dyz", "π")]
    ch = [{"C": a, "M": b, "type": t, "pair_parity": channel(a, b)} for a, b, t in bonds]

    # 3) C–M 跨界面对：C 侧 sp³ 悬挂键 × M 侧 dz²（σ）
    even_w = C_S ** 2 if channel("s", "dz2") == 1 else 0.0
    odd_w = C_PZ ** 2 if channel("pz", "dz2") == -1 else 0.0

    out = {"point_group": "C4v (C(001)/M(001))", "parity_op": "M_z",
           "orbitals": table, "bonds": ch,
           "sp3": {"c_s": C_S, "c_pz": round(C_PZ, 4),
                   "w_s": round(C_S ** 2, 3), "w_pz": round(C_PZ ** 2, 3)},
           "CM_channel": {"sigma_pair_even_weight": round(even_w, 3),
                          "sigma_pair_odd_weight": round(odd_w, 3),
                          "reading": "C 侧 sp³ 悬挂键 × M 侧 dz²：偶分量 25%、奇分量 75%"},
           "criterion": "全 gap（宇称为零）要求偶分量为支配项：w_A1 > 1/2；"
                        "C–M 的自然界面给 1/4，差一半——要 s 主导的界面键或重构"}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1)[:1200])
    print("saved", OUT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_mode.py —— Zd013：轻链-重链的界面局域模

模型：轻链（m, K）n>=0，重链（M, K'）n<=-1，界面弹簧 K_I。
局域模（在轻链带之上）：u_n = u0 alpha^n，v_n = v_-1 beta^{-(n+1)}。
本征方程：P Q = K_I (P + Q)，P = m w^2 - K(1-alpha)，Q = M w^2 - K'(1-beta)，
alpha = (2-A)/2 + sqrt((2-A)^2/4 - 1)，A = m w^2/K（同理 beta, B）。
单位：w0 = sqrt(K/m)，轻链带顶 = 2 w0。
输出：data/interface_mode.json
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_mode.json")


def alpha_of(A):
    """|alpha|<1 的根；A>4 时返回负实数，A<4 时返回 None（传播模）。"""
    if A <= 4:
        return None
    s = math.sqrt((2 - A) ** 2 / 4 - 1)
    return (2 - A) / 2 + s


def F(w2, m, K, M, Kp, KI):
    """w2 = omega^2 / (K/m) 的归一化；返回 F。"""
    A = m * w2 / K
    B = M * w2 / Kp
    a = alpha_of(A)
    b = alpha_of(B)
    if a is None or b is None:
        return None
    P = m * w2 - K * (1 - a)
    Q = M * w2 - Kp * (1 - b)
    return P * Q - KI * (P + Q), P, Q, a, b


def solve(m=1.0, K=1.0, M=1.0, Kp=1.0, KI=2.0):
    """在轻链带顶之上找根；返回 omega/sqrt(K/m) 与 alpha,beta。"""
    lo = 4.0 + 1e-9
    hi = 200.0
    flo = F(lo, m, K, M, Kp, KI)
    fhi = F(hi, m, K, M, Kp, KI)
    if flo is None or fhi is None:
        return None
    if flo[0] * fhi[0] > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = F(mid, m, K, M, Kp, KI)
        if fm is None:
            lo = mid
            continue
        if flo[0] * fm[0] <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    w = math.sqrt(0.5 * (lo + hi))
    return w


def main():
    rows = []
    for ratio in (1.0, 2.0, 3.0, 5.0, 8.0, 12.0):
        for mass in (1.0, 5.0, 20.0, 50.0, 184.0, 1e6):
            w = solve(M=mass, KI=ratio)
            rows.append({"KI_over_K": ratio, "M_over_m": mass,
                         "w_over_w0": round(w, 3) if w else None,
                         "above_band": bool(w and w > 2.0)})
    # 线性拟合 M->inf 的解析式：w^2 = KI^2/(KI-K)（K=m=1）
    analytic = [{"KI_over_K": r, "w_over_w0": round(math.sqrt(r * r / (r - 1)), 3)}
                for r in (2.0, 3.0, 5.0, 8.0, 12.0)]
    out = {"rows": rows, "analytic_Minf": analytic}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"{'KI/K':>6} " + " ".join(f"{'M/m='+str(m):>10}" for m in (1, 5, 20, 50, 184, 'inf')))
    for r in (1.0, 2.0, 3.0, 5.0, 8.0, 12.0):
        vals = []
        for m in (1.0, 5.0, 20.0, 50.0, 184.0, 1e6):
            w = solve(M=m, KI=r)
            vals.append(f"{w:.3f}" if w else "  --  ")
        print(f"{r:6.1f} " + " ".join(f"{v:>10}" for v in vals))
    print("analytic M->inf:", [a["w_over_w0"] for a in analytic])
    print("saved", OUT)


if __name__ == "__main__":
    main()

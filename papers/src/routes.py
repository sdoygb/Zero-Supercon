#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes.py —— Zd004 的三张表：同位素、二维范霍夫、平带-相刚度

口径：Tc 用 Zd002 (0.1) 的 Allen–Dynes 形式；mu* 用 Zb072 命题 11.13。
输出：data/routes.json
"""
from __future__ import annotations
import json, math, os
import numpy as np

from constants import N_TOT

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "routes.json")

AKF = (3 * math.pi ** 2) ** (1 / 3)
AB = 0.529177
N0 = N_TOT  # 2·C_N（见 constants.py）
Q2 = 18.97
ME = 1822.888


def mu_of(a):
    e2 = 2 * (a / AB) / AKF ** 2
    return N0 * 4 * math.pi * e2 / (Q2 + 4 * math.pi * e2 * N0)


def mustar(a, M_amu):
    L = AKF * math.sqrt(M_amu * ME / 2)
    m = mu_of(a)
    return m / (1 + m * math.log(L))


def ad(lam, w, ms):
    den = lam - ms * (1 + 0.62 * lam)
    if den <= 0:
        return 0.0
    return w / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def isotope(a=3.0):
    rows = []
    for lam in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        h = ad(lam, 1.0, mustar(a, 1.0))
        d = ad(lam, 1.0 / math.sqrt(2), mustar(a, 2.0))
        rows.append({"lam": lam, "ratio": d / h,
                     "alpha": -math.log(d / h) / math.log(2)})
    return rows


def dos2d(sigmas=(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002), N=2048):
    k = 2 * np.pi * np.arange(N) / N
    kx, ky = np.meshgrid(k, k, indexing="ij")
    ev = (-2 * (np.cos(kx) + np.cos(ky))).ravel()
    out = []
    for s in sigmas:
        d = float(np.mean(np.exp(-(ev / s) ** 2 / 2)) / (s * math.sqrt(2 * math.pi)))
        out.append({"sigma": s, "N0": d})
    x = np.log([1 / r["sigma"] for r in out])
    y = np.array([r["N0"] for r in out])
    slope = float(np.polyfit(x, y, 1)[0])
    return out, slope


def flat_band(V=0.5, w=0.1, ms=0.1, n=1.0, mmax=60.0):
    rows = []
    best = None
    for m in np.arange(0.5, mmax + 0.25, 0.5):
        lam = 2 * V * m / math.pi
        tmf = ad(lam, w, ms)
        tbkt = math.pi / 2 * n / (4 * m)
        t = min(tmf, tbkt)
        rows.append({"m": float(m), "lam": lam, "TcMF": tmf, "TBKT": tbkt, "Tc": t})
        if best is None or t > best["Tc"]:
            best = dict(rows[-1])
    return rows, best


def main():
    iso = isotope()
    dos, slope = dos2d()
    fb, best = flat_band()
    out = {"isotope": iso, "dos2d": dos, "dos2d_slope": slope,
           "flat": best, "flat_scan": fb}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)

    print("isotope (a=3, H->D):")
    for r in iso:
        print(f"  lam={r['lam']:4.1f}  T_D/T_H={r['ratio']:.4f}  alpha={r['alpha']:.3f}")
    print("2D DOS: slope vs ln(1/sigma) =", round(slope, 4), " 1/(2pi^2)=", round(1 / (2 * math.pi ** 2), 4))
    for r in dos:
        print(f"  sigma={r['sigma']:.3f}  N0={r['N0']:.4f}")
    print("flat band optimum:", {k: round(v, 4) for k, v in best.items()})
    for r in fb:
        if abs(r["m"] - 5) < 1e-9 or abs(r["m"] - 10) < 1e-9 or abs(r["m"] - 15) < 1e-9 \
                or abs(r["m"] - 20) < 1e-9 or abs(r["m"] - 30) < 1e-9 or abs(r["m"] - 50) < 1e-9:
            print(f"  m*={r['m']:4.1f} lam={r['lam']:5.2f} TcMF={r['TcMF']:.4f} TBKT={r['TBKT']:.4f} Tc={r['Tc']:.4f}")
    print("saved", OUT)


if __name__ == "__main__":
    main()

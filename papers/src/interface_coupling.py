#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""interface_coupling.py —— Zd018：界面 g 的单模标定与多模增强

g = D_I q u0, u0 = sqrt(hbar/(2 m omega)), q = 1/a
g_eff = sqrt(M) g; g_req ∝ N_B^{-1/2}
输出：data/interface_coupling.json
"""
from __future__ import annotations
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "interface_coupling.json")
HBAR = 1.0546e-34
AMU = 1.66054e-27
M_C = 12 * AMU
W_EV = 0.15
A_M = 3e-10


def main():
    w = W_EV / 6.582119569e-16
    u0 = math.sqrt(HBAR / (2 * M_C * w))
    q = 1 / A_M
    rows = []
    for D in (2.0, 4.5, 10.0):
        g = D * q * u0
        rows.append({"D_I_eV": D, "g_eV": round(g, 4),
                     "ratio_to_0.43": round(0.43 / g, 1),
                     "ratio_to_0.63": round(0.63 / g, 1)})
    gmin, gmax = rows[0]["g_eV"], rows[-1]["g_eV"]
    meff = [{"M": M, "g_eff_min": round(gmin * math.sqrt(M), 3),
             "g_eff_max": round(gmax * math.sqrt(M), 3)} for M in (10, 30, 100)]
    nb = [{"N_B": n, "scale": round(1 / math.sqrt(n), 2)} for n in (0.2, 0.5, 1.0)]
    out = {"u0_A": round(u0 * 1e10, 4), "q_u0": round(q * u0, 5),
           "rows": rows, "multi": meff, "NB": nb}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"u0 = {out['u0_A']} A, q*u0 = {out['q_u0']}")
    for r in rows:
        print(f"  D={r['D_I_eV']:5.1f} eV  g={r['g_eV']:.4f} eV  ratio {r['ratio_to_0.43']}x / {r['ratio_to_0.63']}x")
    for m in meff:
        print(f"  M={m['M']:4d}  g_eff={m['g_eff_min']}-{m['g_eff_max']} eV")
    for n in nb:
        print(f"  N_B={n['N_B']}  scale={n['scale']}")
    print("saved", OUT)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Zd018 §5：界面 slab 的模型 alpha^2F（4x4x16，8 层 C + 8 层 W，中心弹簧）
# ---------------------------------------------------------------------------
def slab_tc(kappa, DI=5.0, NB=1.0, wcut=300.0, NC=8, NW=8, Nx=4, Ny=4, omega_C=None, m_A=12.0):
    import numpy as np
    Nz = NC + NW
    N = Nx * Ny * Nz
    mCu, mWu = float(m_A), 184.0
    mC = mCu * 1.66054e-27
    idx, atoms = {}, []
    for iz in range(Nz):
        for iy in range(Ny):
            for ix in range(Nx):
                idx[(ix, iy, iz)] = len(atoms); atoms.append((ix, iy, iz))
    mass = np.array([mCu if a[2] < NC else mWu for a in atoms])
    D = np.zeros((3 * N, 3 * N))
    for (ix, iy, iz), i in idx.items():
        for d, (dx, dy, dz) in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            jx, jy, jz = ix + dx, iy + dy, iz + dz
            if dz == 1 and jz >= Nz:
                continue
            jx %= Nx; jy %= Ny
            j = idx[(jx, jy, jz)]
            K = 1.0 if (iz < NC and jz < NC) else (1.0 if (iz >= NC and jz >= NC) else kappa)
            e = np.zeros(3); e[d] = 1.0; blk = K * np.outer(e, e)
            D[3*i:3*i+3, 3*i:3*i+3] += blk / mass[i]
            D[3*j:3*j+3, 3*j:3*j+3] += blk / mass[j]
            D[3*i:3*i+3, 3*j:3*j+3] -= blk / math.sqrt(mass[i] * mass[j])
            D[3*j:3*j+3, 3*i:3*i+3] -= blk / math.sqrt(mass[i] * mass[j])
    w, v = np.linalg.eigh(D)
    lamC = 4.0 / mCu
    if omega_C is None:
        omega_C = 1336.0 * math.sqrt(12.0 / mCu)
    iface = [idx[(ix, iy, NC - 1)] for ix in range(Nx) for iy in range(Ny)]
    kB = 1.380649e-23; hbar = 1.0546e-34; q = 1 / 3e-10
    wK = omega_C * np.sqrt(np.clip(w, 0, None) / lamC)
    S = SL = 0.0
    for i in range(len(w)):
        if wK[i] <= wcut:
            continue
        W = sum(v[3*j:3*j+3, i] @ v[3*j:3*j+3, i] for j in iface)
        om = wK[i] * kB / hbar
        u = math.sqrt(hbar / (2 * mC * om)) * math.sqrt(W)
        g = DI * q * u
        A = g * g / (wK[i] * 8.617333e-5)
        S += A; SL += A * math.log(wK[i])
    lam = 2 * NB * S
    wlog = math.exp(SL / S)
    ms = 0.12
    den = lam - ms * (1 + 0.62 * lam)
    tc = wlog / 1.2 * math.exp(-1.04 * (1 + lam) / den) if den > 0 else 0.0
    return {"kappa": kappa, "lam": lam, "wlog": wlog, "Tc": tc}

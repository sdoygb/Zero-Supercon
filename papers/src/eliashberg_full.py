#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eliashberg_full.py —— Zd022：全谱线性化 Eliashberg 与两矩公式的差

问题：$T_c$ 的两矩公式（Allen–Dynes，用 $\\lambda$ 与 $\\omega_{\\log}$）只用了谱的两个数；
真实 $\\alpha^2F$ 的谱形（峰位、双峰、高频尾）进不去。本脚本用全谱解线性化 Eliashberg，
与同一 $(\\lambda,\\omega_{\\log},\\mu^*)$ 的两矩公式比，量出"谱形缺口"。

方程（各向同性、线性化、$T_c$ 处 $\\Delta\\to0$）：
    Z_n = 1 + (πT/ω_n) Σ_m λ(m-n) sgn(ω_m)
    Δ_n Z_n = πT Σ_m [λ(m-n)-μ*] Δ_m/|ω_m|
    λ(m-n) = ∫ 2ω α²F(ω) / (ω² + (ω_m-ω_n)²) dω
$T_c$ 取"最大本征值 = 1"。矩阵写成 A = diag(πT/Z)·(L-μ*)·diag(1/|ω|)，λ_max(A)=1。

基线：同一 $(\\lambda,\\omega_{\\log})$ 的 Debye 谱（α²F∝ω²）也过一遍同一求解器，
故"材料比值 / Debye 基线"才是净谱形效应（两矩公式本身的拟合偏差被基线吸收）。

数据：Zero 库 `_staging/材料库/jarvis_epc_data_figshare_1058.json`（1058 条 α²F）。
输出：data/eliashberg_full.json
"""
from __future__ import annotations
import json
import math
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "eliashberg_full.json")
SRC = "/Users/oygb/Downloads/Zero/_staging/材料库/jarvis_epc_data_figshare_1058.json"

KB = 0.0861733          # meV/K
MEV_K = 1.0 / KB        # K/meV
MU = 0.10
N_MAT = 150             # Matsubara 半宽
T_MIN, T_MAX = 0.5, 3000.0   # 低于 0.5 K 记为下界（Matsubara 截断失效）
N_BIS = 24


def moments(x, y):
    lam = 2 * np.trapezoid(y / x, x)
    if lam <= 1e-6:
        return lam, float("nan")
    wlog = math.exp(np.trapezoid(2 * y / x * np.log(x), x) / lam)
    return lam, wlog


def debye_spectrum(lam, wlog, n=600):
    """同一 (λ, ω_log) 的 Debye 谱：α²F=Cω², ω<ω_D；ω_log=ω_D e^{-1/2}。"""
    wD = wlog * math.exp(0.5)
    x = np.linspace(1e-3, wD, n)
    y = lam / wD ** 2 * x ** 2
    return x, y


def kernel_table(step, nmax, x, y):
    ks = np.arange(-nmax, nmax + 1)
    d = ks * step
    W = (2 * x * y)[None, :] / (x[None, :] ** 2 + d[:, None] ** 2)
    return np.trapezoid(W, x, axis=1)


def tc_full(x, y, mu=MU, N=N_MAT, tmin=None, tmax=None):
    wcmu = 10.0 * x.max()

    def ev(T):
        Tm = T * KB
        n = np.arange(-N, N + 1)
        wn = (2 * n + 1) * np.pi * Tm
        step = wn[1] - wn[0]
        g = kernel_table(step, 2 * N, x, y)
        L = g[(n[:, None] - n[None, :]) + 2 * N]
        mask = (np.abs(wn) < wcmu).astype(float)
        Z = 1 + (np.pi * Tm / wn) * np.sum(L * np.sign(wn)[None, :], axis=1)
        if np.any(Z <= 0.2):
            return np.nan   # Matsubara 截断未收敛
        A = np.diag(np.pi * Tm / Z) @ (L - mu * np.outer(mask, mask)) @ np.diag(1 / np.abs(wn))
        return float(np.linalg.eigvals(A).real.max())

    t_lo = T_MIN if tmin is None else tmin
    t_hi = T_MAX if tmax is None else tmax
    f1, f2 = ev(t_lo), ev(t_hi)
    if not np.isfinite(f1) or not np.isfinite(f2):
        return float("nan")
    if f1 < 1:
        return 0.0
    if f2 > 1:
        return t_hi
    lo, hi = math.log(t_lo), math.log(t_hi)
    for _ in range(N_BIS):
        mid = 0.5 * (lo + hi)
        f = ev(math.exp(mid))
        if not np.isfinite(f):
            return float("nan")
        if f > 1:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def ad(lam, wlog_meV, mu=MU):
    den = lam - mu * (1 + 0.62 * lam)
    if den <= 1e-9:
        return float("nan")
    return wlog_meV / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def shape(x, y, lam):
    """谱形描述：峰位、峰位/ω_log、按 λ 加权的高频份额（ω>ω_log）。"""
    wlog = math.exp(np.trapezoid(2 * y / x * np.log(x), x) / lam)
    wpeak = float(x[int(np.argmax(y))])
    w = 2 * y / x
    hi = x > wlog
    frac_hi = float(np.trapezoid(w[hi], x[hi]) / np.trapezoid(w, x))
    return {"wpeak_meV": round(wpeak, 3), "wpeak_over_wlog": round(wpeak / wlog, 3),
            "frac_lambda_above_wlog": round(frac_hi, 3)}


BASE_LAMS = np.geomspace(0.15, 4.0, 20)
_BASE = None


def base_of(lam):
    """Debye 基线 tE/tA，只依赖 λ（两式都随 ω_log 线性缩放）。"""
    global _BASE
    if _BASE is None:
        _BASE = []
        for l0 in BASE_LAMS:
            xd, yd = debye_spectrum(l0, 20.0)
            ld, wd = moments(xd, yd)
            v = ad(ld, wd) * MEV_K
            _BASE.append(tc_full(xd, yd) / v if v > 0 else float('nan'))
        _BASE = np.array(_BASE, float)
        ok = np.isfinite(_BASE)
        _BASE = np.interp(np.arange(len(_BASE)), np.arange(len(_BASE))[ok], _BASE[ok])
    return float(np.interp(math.log(lam), np.log(BASE_LAMS), _BASE))




def bcs_limit(lam=0.2, wc=30.0, mu=0.0):
    """自检：常数核 + 尖锐 cutoff 的 BCS 极限（应给 ~1.134 wc e^{-1/lam}）。"""
    def ev(T):
        Tm = T * KB
        N = int(np.ceil(wc / (2 * np.pi * Tm))) + 2
        n = np.arange(-N, N + 1)
        wn = (2 * n + 1) * np.pi * Tm
        mask = (np.abs(wn) <= wc).astype(float)
        L = lam * np.outer(mask, mask)
        Z = 1 + (np.pi * Tm / wn) * np.sum(L * np.sign(wn)[None, :], axis=1)
        A = np.diag(np.pi * Tm / Z) @ (L - mu * np.outer(mask, mask)) @ np.diag(1 / np.abs(wn))
        return float(np.linalg.eigvals(A).real.max())

    lo, hi = math.log(0.5), math.log(20000.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if ev(math.exp(mid)) > 1:
            lo = mid
        else:
            hi = mid
    T = math.exp(0.5 * (lo + hi))
    ref = 1.134 * (wc / KB) * math.exp(-1 / lam)
    return {"Tc_K": round(T, 2), "ref_BCS_K": round(ref, 2), "ratio": round(T / ref, 3)}


def main():
    d = json.load(open(SRC))
    rows = []
    skipped = {"bad_lambda": 0, "bad_wlog": 0, "moments_mismatch": 0, "solver_nan": 0}
    for r in d:
        try:
            x = np.asarray(r["a2F_original_x"], float)
            y = np.asarray(r["a2F_original_y"], float)
        except Exception:
            skipped["bad_lambda"] += 1
            continue
        ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y >= 0)
        x, y = x[ok], y[ok]
        if len(x) < 10:
            skipped["bad_lambda"] += 1
            continue
        lam, wlog = moments(x, y)
        if not np.isfinite(lam) or lam <= 0.02:
            skipped["bad_lambda"] += 1
            continue
        if not np.isfinite(wlog) or not (10 < wlog * MEV_K < 3000):
            skipped["bad_wlog"] += 1
            continue
        jl = r.get("lamb", float("nan"))
        if np.isfinite(jl) and abs(lam - jl) / max(jl, 1e-9) > 0.10:
            skipped["moments_mismatch"] += 1
            continue
        tE = tc_full(x, y)
        if not np.isfinite(tE):
            skipped["solver_nan"] += 1
            continue
        tA = ad(lam, wlog)
        base = base_of(lam)
        rows.append({"jid": r["jid"], "lam": round(lam, 4), "wlog_K": round(wlog * MEV_K, 1),
                     "Tc_full": round(tE, 3), "Tc_AD": round(tA * MEV_K, 3),
                     "ratio": round(tE / (tA * MEV_K), 3) if tA > 0 else None,
                     "debye_base": round(base, 3) if np.isfinite(base) else None,
                     "shape_effect": round((tE / (tA * MEV_K)) / base, 3) if (tA > 0 and np.isfinite(base) and base > 0) else None,
                     **shape(x, y, lam), "stability": r.get("stability")})
    # 统计
    def stats(vs):
        vs = [v for v in vs if v is not None and np.isfinite(v)]
        if not vs:
            return {}
        vs = sorted(vs)
        q = lambda p: vs[min(len(vs) - 1, int(p * len(vs)))]
        return {"n": len(vs), "median": round(q(0.5), 3), "q25": round(q(0.25), 3),
                "q75": round(q(0.75), 3), "min": round(vs[0], 3), "max": round(vs[-1], 3)}

    bins = [(0.02, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 1.5), (1.5, 5.0)]
    by_bin = {}
    for a, b in bins:
        sub = [r for r in rows if a <= r["lam"] < b]
        by_bin[f"{a}-{b}"] = {"n": len(sub),
                              "ratio": stats([r["ratio"] for r in sub]),
                              "shape_effect": stats([r["shape_effect"] for r in sub])}
    out = {"n_total": len(d), "n_used": len(rows), "skipped": skipped,
           "mu_star": MU, "N_matsubara": N_MAT,
           "ratio_stats": stats([r["ratio"] for r in rows]),
           "shape_effect_stats": stats([r["shape_effect"] for r in rows]),
           "by_bin": by_bin,
           "rows": rows}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("used", len(rows), "skipped", skipped)
    print("ratio", out["ratio_stats"])
    print("shape_effect", out["shape_effect_stats"])
    for k, v in by_bin.items():
        print(k, v["n"], v["shape_effect"])
    print("saved", OUT)


if __name__ == "__main__":
    main()

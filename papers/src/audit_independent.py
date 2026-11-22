#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_independent.py —— Zd001–Zd019 推导与计算的独立复核（不复用 src/ 里的实现）

口径：所有公式在本文件里重写一遍；只为与正文/JSON 对数，不改变任何正文。
输出：data/audit_independent.json + 屏幕摘要
"""
import csv, json, math, statistics as st, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CJ = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"
CM = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"
KBEV = 8.617333262e-5


def AD(lam, w, ms):
    den = lam - ms * (1 + 0.62 * lam)
    return 0.0 if den <= 0 else w / 1.2 * math.exp(-1.04 * (1 + lam) / den)


def lam_req(w, ms=0.12, T=300.0):
    rhs = math.log(w / (1.2 * T))
    return (1.04 + rhs * ms) / (rhs * (1 - 0.62 * ms) - 1.04)


def ceiling(ms):
    return (1 / 1.2) * math.exp(-1.04 / (1 - 0.62 * ms))


out = {}

# ---------- 1. 窗口与天花板 ----------
out["ceiling"] = {str(ms): round(ceiling(ms), 4) for ms in (0.085, 0.10, 0.12)}
out["w_req_vs_lam"] = {str(ms): {str(l): round(lam_req(1.0, ms) if False else 0, 2)
                                 for l in (2, 3, 5)} for ms in (0.12,)}
# λ_req 的反函数：给定 λ，求 ω
w_req_tab = {}
for ms in (0.085, 0.12):
    row = {}
    for lam in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
        den = lam - ms * (1 + 0.62 * lam)
        row[str(lam)] = round(1.2 * 300 * math.exp(1.04 * (1 + lam) / den), 1) if den > 0 else None
    w_req_tab[str(ms)] = row
out["w_req"] = w_req_tab

# ---------- 2. μ* ----------
AB = 0.529177
AKF = (3 * math.pi ** 2) ** (1 / 3)
N0 = 0.2854  # 2·C_N（见 constants.py）
Q2 = 18.97


def mu_star(a_A, M_amu, q2=Q2, n0=N0):
    e2 = 2 * (a_A / AB) / AKF ** 2
    mu = n0 * 4 * math.pi * e2 / (q2 + 4 * math.pi * e2 * n0)
    L = AKF * math.sqrt((M_amu / 1822.888) / 2.0)
    return mu / (1 + mu * math.log(L)), mu, L


out["mu_star_check"] = {f"a={a},M={M}": [round(x, 4) for x in mu_star(a, M)]
                        for a, M in ((3.0, 1.0), (3.0, 12.0), (2.5, 1.0), (2.0, 195.0))}

# ---------- 3. JARVIS 统计 ----------
rows = list(csv.DictReader(open(CJ)))
for r in rows:
    for k in ("lam", "wlog_K", "a_A", "M_amu", "Tc_K", "lam_a2f", "wlog_a2f_K"):
        try:
            r[k] = float(r[k])
        except Exception:
            r[k] = None
filt = [r for r in rows if r["lam"] > 0.02 and r["wlog_K"] > 1.0]
out["jarvis"] = {
    "n_total": len(rows),
    "n_effective": len(filt),
    "max_wlog": round(max(r["wlog_K"] for r in rows), 1),
    "max_over_1091": sum(1 for r in rows if r["wlog_K"] > 1091),
    "duplicates_same_formula": len(rows) - len(set(r["formula"] for r in rows)),
    "median_lam": round(st.median([r["lam"] for r in filt]), 3),
    "q25_lam": round(st.quantiles([r["lam"] for r in filt], n=4)[0], 3),
    "q75_lam": round(st.quantiles([r["lam"] for r in filt], n=4)[2], 3),
    "frac_gt1": round(sum(1 for r in filt if r["lam"] > 1) / len(filt), 3),
    "frac_gt2": round(sum(1 for r in filt if r["lam"] > 2) / len(filt), 3),
    "max_lam": round(max(r["lam"] for r in filt), 2),
    "sub_lam_gt2_max_wlog": round(max((r["wlog_K"] for r in filt if r["lam"] > 2), default=0), 1),
    "sub_lam_gt5_max_wlog": round(max((r["wlog_K"] for r in filt if r["lam"] > 5), default=0), 1),
    "n_lam_gt2": sum(1 for r in filt if r["lam"] > 2),
    "n_lam_gt5": sum(1 for r in filt if r["lam"] > 5),
}
# 斜率（独立最小二乘）
xs = [math.log(r["lam"]) for r in filt]
ys = [math.log(r["wlog_K"]) for r in filt]
mx, my = st.mean(xs), st.mean(ys)
sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
sxx = sum((x - mx) ** 2 for x in xs)
syy = sum((y - my) ** 2 for y in ys)
out["jarvis"]["slope_lnw_lnlam"] = round(sxy / sxx, 3)
out["jarvis"]["corr"] = round(sxy / math.sqrt(sxx * syy), 3)
out["jarvis"]["alpha_implied"] = round(1 + 1 / (2 * (sxy / sxx)), 3)

# MoN 的倍数
mon = next(r for r in rows if r["formula"] == "MoN")
ms_mon = float(mon["mustar"])
w_need = 1.2 * 300 * math.exp(1.04 * (1 + mon["lam"]) /
                              (mon["lam"] - ms_mon * (1 + 0.62 * mon["lam"])))
out["mon"] = {"lam": mon["lam"], "wlog": mon["wlog_K"], "mustar": ms_mon,
              "factor": round(w_need / mon["wlog_K"], 2),
              "w_needed": round(w_need, 0)}
# α²F 列与 lam/wlog 列的差
d_lam = [abs(r["lam"] - r["lam_a2f"]) for r in rows if r["lam_a2f"]]
d_w = [abs(r["wlog_K"] - r["wlog_a2f_K"]) for r in rows if r["wlog_a2f_K"]]
out["jarvis"]["median_absdiff_lam_vs_a2f"] = round(st.median(d_lam), 3)
out["jarvis"]["median_absdiff_wlog_vs_a2f"] = round(st.median(d_w), 1)

# ---------- 4. Zd005 设计方程 × Zd015 RPA ----------
design = []
for wA_K in (1336, 1400, 1700, 1916, 2000, 2500):
    wA = wA_K * KBEV
    for NB in (0.5, 1.0, 2.0):
        lr = lam_req(wA_K, 0.12)
        g = math.sqrt(lr * wA / (2 * NB))
        lt = 2 * NB * g * g / wA
        # RPA: ω~² = ω² - 4 g² N ω
        soft = 1 - 4 * g * g * NB / wA
        design.append({"wA_K": wA_K, "NB": NB, "g_req": round(g, 3),
                       "g_over_w": round(g / wA, 2), "lam_eff": round(lt, 2),
                       "rpa_soft_ratio": round(soft, 2),
                       "lam_max_stable": 0.5,
                       "Tc_at_stable_max": round(AD(0.5, wA_K, 0.12), 0)})
out["design_vs_rpa"] = design

# ---------- 5. Zd019 复核 ----------
def L(T, w):
    return math.log(1.13 * w / T) if 0 < T < 1.13 * w else -1e9


def L12(T, w, d):
    if d <= 1e-9:
        return L(T, w)
    # 简单梯形积分（2000 点），不用 scipy
    n = 2000
    h = w / n
    s = 0.0
    for i in range(n + 1):
        x = i * h
        f = math.tanh(x / (2 * T)) / (x + d / 2)
        s += f * (0.5 if i in (0, n) else 1.0)
    return s * h


def two_band_tc(l11, w1, l22, w2, l12, w12, ms0=0.12, ms12=0.04, delta=0.0):
    # det M = 1, 二分
    def det(T):
        a = (l11 - ms0) * L(T, w1)
        d = (l22 - ms0) * L(T, w2)
        b = (l12 - ms12) * L12(T, w12, delta)
        return (a * d - b * b) - 1.0
    lo, hi = 0.5, min(w1, w2) * 1.1
    if det(lo) * det(hi) > 0:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if det(lo) * det(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


delta_tab = {}
for name, w in (("mgb2", 846.0), ("diamond", 1916.0)):
    delta_tab[name] = [{"frac": f, "Tc": round(two_band_tc(0.30, w, 0.10, w, 0.20, w, delta=f * w), 1)}
                       for f in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)]
out["zd019_delta"] = delta_tab

# 从正文参数直接算 λ_max（δ=0）
lam_mat = [[0.30 - 0.12, 0.20 - 0.04], [0.20 - 0.04, 0.10 - 0.12]]
ev = sorted([ (lam_mat[0][0]+lam_mat[1][1])/2 + math.sqrt(((lam_mat[0][0]-lam_mat[1][1])/2)**2 + lam_mat[0][1]**2),
              (lam_mat[0][0]+lam_mat[1][1])/2 - math.sqrt(((lam_mat[0][0]-lam_mat[1][1])/2)**2 + lam_mat[0][1]**2)])
out["zd019_lam_max_delta0"] = [round(x, 3) for x in ev]

# 质量匹配表：正文 vs 代码
def neff(eta, wD=0.165, Dsp=0.0056, cap=40):
    return min(cap, math.sqrt(wD / (Dsp * eta))) if eta > 0 else cap
mm = []
for eta, text_tc in ((0.01, 298), (0.03, 263), (0.10, 177), (0.20, 135)):
    n = neff(eta)
    lam_no_p = (0.46 - 0.12) + (n - 1) * (0.05 - 0.04)
    lam_p6 = (0.46 - 0.12) + (n - 1) * (0.05 - 0.04) * 0.6
    mm.append({"eta": eta, "Neff": round(n, 1), "Tc_text": text_tc,
               "Tc_no_parity_double_sub": round(AD(lam_no_p, 1916, 0.12), 1),
               "Tc_parity0.6": round(AD(lam_p6, 1916, 0.12), 1),
               "lam_no_p": round(lam_no_p, 3), "lam_p6": round(lam_p6, 3)})
out["zd019_massmatch"] = mm

# μ* 双减：λ 里已减 MS0，AD 里再减一次
out["zd019_double_sub"] = {
    "best_lambda_with_double": 2.68,
    "Tc_double": round(AD(2.68, 1916, 0.12), 1),
    "Tc_single": round(AD(2.68 + 0.12 - 0.12 * (1 + 0.62 * 2.68) * 0 + 0, 1916, 0.12), 1),
    "note": "λ_tot=(λ0-μ0)+(N-1)net 已含 μ0；AD 又按 μ*(1+0.62λ) 减一次",
}
# 更准确：无 μ0 的 λ 应为 λ_tot+0.12=2.80
out["zd019_double_sub"]["Tc_single_correct"] = round(AD(2.80, 1916, 0.12), 1)

# ---------- 6. Zd004 平带：独立网格 ----------
flat = []
for v, ms in ((0.5, 0.10), (0.5, 0.12)):
    for mstar in (2, 5, 10, 17.5, 25, 50):
        lam = 2 * v * mstar / math.pi
        Tmf = (v / 1.2) * math.exp(-1.04 * (1 + lam) / (lam - ms * (1 + 0.62 * lam))) if lam > ms * (1 + 0.62 * lam) else 0
        Tbkt = math.pi / (8 * mstar)
        flat.append({"mstar": mstar, "lam": round(lam, 2), "Tmf": round(Tmf, 4),
                     "Tbkt": round(Tbkt, 4), "Tmin": round(min(Tmf, Tbkt), 4)})
out["zd004_flat"] = flat

# ---------- 7. 3D 半满 vH 的 DOS 系数（数值，独立） ----------
def dos_cubic(E, eta, n=300):
    # N(E) = (1/(2π)^3) ∫ dk δ_η(E-ε_k) / t ；ε_k/t=-2(cx+cy+cz)
    s = 0.0
    h = 2 * math.pi / n
    for i in range(n):
        cx = math.cos(-math.pi + (i + 0.5) * h)
        for j in range(n):
            cy = math.cos(-math.pi + (j + 0.5) * h)
            for k in range(n):
                cz = math.cos(-math.pi + (k + 0.5) * h)
                e = -2 * (cx + cy + cz)
                s += (eta / math.pi) / ((e - E) ** 2 + eta ** 2)
    return s * h ** 3 / (2 * math.pi) ** 3


out["vhs_dos"] = {f"eta={eta}": round(dos_cubic(0.0, eta, 120), 4) for eta in (0.02, 0.05, 0.1)}

json.dump(out, open(os.path.join(ROOT, "data", "audit_independent.json"), "w"),
          ensure_ascii=False, indent=1)

print("=== 天花板 ===", out["ceiling"])
print("=== ω_req ===", json.dumps(w_req_tab, ensure_ascii=False))
print("=== μ* ===", json.dumps(out["mu_star_check"], ensure_ascii=False))
print("=== JARVIS ===", json.dumps(out["jarvis"], ensure_ascii=False))
print("=== MoN ===", out["mon"])
print("=== 设计方程 vs RPA ===")
for d in design:
    print(d)
print("=== Zd019 δ 表 ===", json.dumps(delta_tab, ensure_ascii=False))
print("=== Zd019 λmax(δ=0) ===", out["zd019_lam_max_delta0"])
print("=== Zd019 质量匹配 ===")
for m in mm:
    print(m)
print("=== Zd019 μ* 双减 ===", out["zd019_double_sub"])
print("=== Zd004 平带 ===")
for f in flat:
    print(f)
print("=== vH DOS ===", out["vhs_dos"])

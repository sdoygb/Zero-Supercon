#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_claims.py —— 把 Zd 系列正文里的头条数字逐条重算一遍

每行：claim 名 / 正文值 / 重算值 / 判（|差| <= 容差）
输出：data/verify_claims.json（并打印 PASS/FAIL 表）
"""
from __future__ import annotations
import csv, json, math, os
import numpy as np

from constants import C_LAM_G1, C_LAM_G13, C_N, N_TOT

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JARVIS = "/Users/oygb/Downloads/Zero/_staging/批量Tc_全表.csv"
MP = "/Users/oygb/Downloads/Zero/_staging/MP_声子_分档_全表.csv"
OUT = os.path.join(ROOT, "data", "verify_claims.json")

AKF = (3 * math.pi ** 2) ** (1 / 3)
AB = 0.529177
N0 = N_TOT  # 2·C_N（见 constants.py）
Q2 = 18.97
ME = 1822.888
HBAR = 6.582119569e-16
KB = 8.617333e-5


def ceiling(ms):
    return (1 / 1.2) * math.exp(-1.04 / (1 - 0.62 * ms))


def wreq(lam, ms=0.1, T=300.0):
    den = lam - ms * (1 + 0.62 * lam)
    return 1.2 * T * math.exp(1.04 * (1 + lam) / den)


def lam_req(w, ms=0.1, T=300.0):
    rhs = math.log(w * (1 / 1.2) / T)
    return (1.04 + rhs * ms) / (rhs * (1 - 0.62 * ms) - 1.04)


def rows(path):
    out = []
    for r in csv.DictReader(open(path)):
        try:
            out.append({"w": float(r["wlog_K"]), "l": float(r["lam"]),
                        "tc": float(r["Tc_K"]) if "Tc_K" in r and r["Tc_K"] else None,
                        "M": float(r["M_amu"]) if "M_amu" in r else None,
                        "a": float(r["a_A"]) if "a_A" in r else None,
                        "f": r.get("formula", "")})
        except Exception:
            pass
    return out


def _st_median(xs):
    import statistics
    return statistics.median(xs)


def main():
    j = [r for r in rows(JARVIS) if r["w"] > 1 and r["l"] > 0.02]
    m = [r for r in rows(MP) if r["w"] > 1 and r["l"] > 0.02]
    claims = []

    def ck(name, text, calc, tol=0.01, abs_tol=0.0):
        ok = abs(text - calc) <= max(tol * max(1.0, abs(calc)), abs_tol)
        claims.append({"claim": name, "text": text, "calc": round(calc, 4),
                       "tol": tol, "pass": bool(ok)})

    ck("ceiling mu*=0.077", 0.280, ceiling(0.077), 0.002)
    ck("ceiling mu*=0.10", 0.275, ceiling(0.10), 0.002)
    ck("ceiling mu*=0.113", 0.272, ceiling(0.113), 0.002)
    ck("w_min mu*=0.077", 1073, 300 / ceiling(0.077), 0.002)
    ck("w_min mu*=0.10", 1091, 300 / ceiling(0.10), 0.002)
    ck("w_min mu*=0.113", 1101, 300 / ceiling(0.113), 0.002)
    ck("c_N per spin", 0.1427, C_N, 0.005)
    ck("w_req lam=2", 2086, wreq(2.0), 0.002)
    ck("w_req lam=5", 1402, wreq(5.0), 0.002)
    ck("w_req inf", 1091, 300 / ceiling(0.10), 0.002)
    ck("JARVIS max wlog", 912.4, max(r["w"] for r in j), 0.001)
    ck("JARVIS max wlog lam>=1", 307.3, max(r["w"] for r in j if r["l"] >= 1), 0.001)
    ck("JARVIS max wlog lam>=2", 172.6, max(r["w"] for r in j if r["l"] >= 2), 0.001)
    ck("JARVIS max Tc", 33.0, max(r["tc"] for r in j if r["tc"] is not None), 0.01)
    # corr
    x = np.log([r["l"] for r in j]); y = np.log([r["w"] for r in j])
    ck("JARVIS corr", -0.775, float(np.corrcoef(x, y)[0, 1]), 0.01)
    ck("JARVIS slope", -0.806, float(np.polyfit(x, y, 1)[0]), 0.01)
    x2 = np.log([r["l"] for r in m]); y2 = np.log([r["w"] for r in m])
    ck("MP slope", -0.884, float(np.polyfit(x2, y2, 1)[0]), 0.01)
    ck("MP n", 7678, len(m), 0.001)
    wmin = 300 / ceiling(0.1)
    hi = [r for r in m if r["w"] > wmin]
    ck("MP n above wmin", 15, len(hi), 0.001)
    ck("MP max lam above wmin", 0.270, max(r["l"] for r in hi), 0.01)
    inside = [r for r in hi if r["l"] >= lam_req(r["w"])]
    ck("MP n inside", 0, len(inside), 0.001)
    closest = max(hi, key=lambda r: r["l"] / lam_req(r["w"]))
    ck("MP closest ratio (max)", 0.042, closest["l"] / lam_req(closest["w"]), 0.01, 0.002)

    # isotope
    def ms(a, M):
        L = AKF * math.sqrt(M * ME / 2)
        mu = N0 * 4 * math.pi * (2 * (a / AB) / AKF ** 2) / (Q2 + 4 * math.pi * (2 * (a / AB) / AKF ** 2) * N0)
        return mu / (1 + mu * math.log(L))

    def ad(lam, w, mstar):
        den = lam - mstar * (1 + 0.62 * lam)
        return w / 1.2 * math.exp(-1.04 * (1 + lam) / den)

    for lam, want in ((0.5, 0.428), (5.0, 0.494)):
        h = ad(lam, 1, ms(3, 1)); d = ad(lam, 1 / math.sqrt(2), ms(3, 2))
        ck(f"isotope alpha lam={lam}", want, -math.log(d / h) / math.log(2), 0.01)

    # layer route
    ck("g_req 1700K N=1 (mu*=0.12)", 0.476, math.sqrt(lam_req(1700, ms=0.12) * (1700 * 8.617333e-5) / 2), 0.01)
    # Zd009 natural coupling
    import statistics as _st
    C = C_N * (3 * math.pi ** 2) ** (2 / 3)
    ck("C = cN kappa^2", 1.37, C, 0.01)
    lams = sorted(r["l"] for r in j)
    ck("JARVIS median lam", 0.465, _st.median(lams), 0.01)
    ck("JARVIS frac lam>1", 0.258, sum(1 for x in lams if x > 1) / len(lams), 0.02)
    ck("JARVIS frac lam>2", 0.124, sum(1 for x in lams if x > 2) / len(lams), 0.02)
    ck("JARVIS frac lam>5", 0.029, sum(1 for x in lams if x > 5) / len(lams), 0.02)
    ck("D/t at lam=0.465", 0.58, math.sqrt(0.465 / C), 0.02)
    ck("D/t at lam=2", 1.21, math.sqrt(2 / C), 0.01)

    # Zd010 candidate scan
    ck("JARVIS n above wmin", 0, sum(1 for r in j if r["w"] > 1091), 0.001)
    ck("JARVIS lam>5 max wlog", 48.7, max(r["w"] for r in j if r["l"] > 5), 0.01)
    ck("lam0 B2Mg", 5.54, 0.696 * (846.5 / 300) ** 2, 0.01)
    ck("lam0 B2CN", 5.14, 0.555 * (912.4 / 300) ** 2, 0.01)

    # Zd011 pair table
    def greq(w, n=1.0, ms=0.12):
        return math.sqrt(lam_req(w, ms=ms) * (w * 8.617333e-5) / (2 * n))
    ck("g_req BiCCaO N=1", 0.43, greq(2240.9), 0.01)
    ck("g_req C 1336 N=1", 0.63, greq(1336.5), 0.01)
    ck("g_req CHN N=1", 0.59, greq(1377.5), 0.01)
    ck("lam0 MoN", 1.43, 1.367 * (307.3 / 300) ** 2, 0.02)
    ck("lam0 V", 0.72, 0.959 * (260.2 / 300) ** 2, 0.02)

    # Zd012 interface phonons
    ck("CTi wlog", 713.0, max(r["w"] for r in j if r["f"] == "CTi"), 0.01)
    ck("MoN wlog (JARVIS)", 307.3, max(r["w"] for r in j if r["f"] == "MoN"), 0.01)
    ck("diamond wlog", 1336.5, max(r["w"] for r in m if r["f"] == "C"), 0.01)

    # Zd013 interface localized mode
    import importlib.util as _iu
    _sp = _iu.spec_from_file_location("im", os.path.join(ROOT, "src", "interface_mode.py"))
    _im = _iu.module_from_spec(_sp); _sp.loader.exec_module(_im)
    ck("w_I (KI=3, M=184)", 2.123, _im.solve(M=184.0, KI=3.0), 0.01)
    ck("w_I (KI=5, M=inf)", 2.500, _im.solve(M=1e6, KI=5.0), 0.01)
    ck("w_I (KI=2, M=184)", 2.000, _im.solve(M=184.0, KI=2.0), 0.01)

    # Zd014 anharmonicity
    _CAP = C_N * (3 * math.pi ** 2) ** (2 / 3)
    ck("C (gL=1/3)", 0.455, _CAP / 3, 0.01)
    ck("C (gL=1)", 1.37, _CAP, 0.01)
    ck("n_eff median (gL=1/3)", 1.009, _st.median([math.sqrt(l / (_CAP / 3)) for l in lams]), 0.01)
    ck("n at lam=2 (gL=1/3)", 2.10, math.sqrt(2 / (_CAP / 3)), 0.01)

    # Zd015 stability
    import json as _json
    _d = _json.load(open("/Users/oygb/Downloads/Zero/_staging/材料库/jarvis_epc_data_figshare_1058.json"))
    _st = [x["lamb"] for x in _d if x.get("stability") == "stable" and x["lamb"] > 0]
    _un = [x["lamb"] for x in _d if x.get("stability") == "unstable" and x["lamb"] > 0]
    ck("stable median lam", 0.359, _st_median(_st), 0.01)
    ck("unstable median lam", 1.107, _st_median(_un), 0.01)
    ck("stable frac lam>2", 0.018, sum(1 for x in _st if x > 2) / len(_st), 0.02)
    ck("joint lam>2 w>1091 (stable)", 0,
       sum(1 for x in _d if x.get("stability") == "stable" and x["lamb"] > 2 and x["wlog"] > 1091), 0.001)

    # Zd017 interface lifetime
    tau0 = 2 * math.pi * HBAR / 0.2
    ck("tau(1.5 eV)/1e11", 3.27, tau0 * math.exp(1.5 / (KB * 300)) / 1e11, 0.01)
    ck("tau(2.0 eV)/1e19", 8.20, tau0 * math.exp(2.0 / (KB * 300)) / 1e19, 0.01)

    # Zd018 interface coupling
    _u0 = math.sqrt(1.0546e-34 / (2 * 12 * 1.66054e-27 * (0.15 / 6.582119569e-16)))
    ck("u0 (A)", 0.0341, _u0 * 1e10, 0.01)
    ck("g(D=10 eV)", 0.1136, 10 * (1 / 3e-10) * _u0, 0.01)
    ck("g_eff(M=30, gmax)", 0.622, 30 ** 0.5 * 10 * (1 / 3e-10) * _u0, 0.01)

    # Zd018 §5 slab model
    _sp2 = _iu.spec_from_file_location("ic", os.path.join(ROOT, "src", "interface_coupling.py"))
    _ic = _iu.module_from_spec(_sp2); _sp2.loader.exec_module(_ic)
    _sl = _ic.slab_tc(3.0)
    ck("slab lam(kappa=3)", 4.14, float(_sl["lam"]), 0.01)
    ck("slab Tc(kappa=3)", 208.0, float(_sl["Tc"]), 0.01)
    _s2 = _ic.slab_tc(5.0, DI=15.0, omega_C=1500.0)
    ck("slab Tc(wC=1500,DI=15)", 303.0, float(_s2["Tc"]), 0.01)
    _s3 = _ic.slab_tc(5.0, DI=10.0, omega_C=2000.0)
    ck("slab Tc(wC=2000,DI=10)", 320.0, float(_s3["Tc"]), 0.01)
    _s4 = _ic.slab_tc(5.0, DI=15.0, m_A=6.9)
    ck("slab Tc(mA=6.9,DI=15)", 371.0, float(_s4["Tc"]), 0.01)
    _s5 = _ic.slab_tc(5.0, DI=15.0, m_A=1.0)
    ck("slab Tc(mA=1,DI=15)", 954.0, float(_s5["Tc"]), 0.01)

    # Zd018 §6 real-frequency optical model
    _hbar = 6.582119569e-16
    _mC = 12 * 1.66054e-27
    _woptK = 1916.0
    _wopt = _woptK * 8.617333e-5
    _uopt = math.sqrt(1.0546e-34 / (2 * _mC * (_wopt / _hbar)))
    def _lam_opt(DI, NM):
        g = DI * (1 / 3e-10) * _uopt
        return 2 * NM * g * g / (_woptK * 8.617333e-5)
    def _tc_opt(lam, w=1916.0, ms=0.12):
        den = lam - ms * (1 + 0.62 * lam)
        return w / 1.2 * math.exp(-1.04 * (1 + lam) / den) if den > 0 else 0.0
    ck("lam_opt DI=10 NM=30", 4.26, _lam_opt(10, 30), 0.01)
    ck("Tc_opt DI=10 NM=30", 382.0, _tc_opt(_lam_opt(10, 30)), 0.01)
    _cross = next(n for n in range(1, 200) if _tc_opt(_lam_opt(10, n)) >= 300)
    ck("cross NM (DI=10)", 18, _cross, 0.001)

    # ---- Zd019: 两带能隙矩阵与设计窗口 ----
    import numpy as _np
    from scipy.integrate import quad as _quad
    def _L(T, w):
        return math.log(1.13 * w / T) if 0 < T < 1.13 * w else -1e9
    def _L12(T, w, d):
        if d <= 0:
            return _L(T, w)
        return _quad(lambda x: math.tanh(x / (2 * T)) / (x + d / 2), 0, w, limit=200)[0]
    def _tb(l11, w1, l22, w2, l12, w12, ms0=0.12, ms12=0.04, d=0.0):
        def eig(T):
            M = _np.array([[(l11 - ms0) * _L(T, w1), (l12 - ms12) * _L12(T, w12, d)],
                           [(l12 - ms12) * _L12(T, w12, d), (l22 - ms0) * _L(T, w2)]])
            return _np.linalg.eigvalsh(M)[-1]
        Ts = _np.logspace(math.log10(0.2), math.log10(max(w1, w2) * 1.1), 800)
        roots = [Ts[i] for i in range(1, len(Ts)) if (eig(Ts[i]) - 1) * (eig(Ts[i - 1]) - 1) < 0]
        return max(roots) if roots else 0.0
    def _ad(lam, w=1916.0, ms=0.12):
        den = lam - ms * (1 + 0.62 * lam)
        return w / 1.2 * math.exp(-1.04 * (1 + lam) / den) if den > 0 else 0.0
    ck("lam_req(1916)", 2.445, lam_req(1916.0, ms=0.12), 0.002)
    ck("lam_req(5180)", 0.952, lam_req(5180.0, ms=0.12), 0.002)
    ck("Zd019 two-band delta=0 diamond", 53.0, _tb(0.3, 1916.0, 0.1, 1916.0, 0.2, 1916.0), 0.03)
    ck("Zd019 two-band delta=383 diamond", 17.0, _tb(0.3, 1916.0, 0.1, 1916.0, 0.2, 1916.0, d=383.2), 0.05)
    ck("Zd019 mu*=0.04 lam12=0.05", 63.0, _ad((0.46 - 0.12) + 39 * (0.05 - 0.04)), 0.05)
    ck("Zd019 best scenario", 316.0, _ad((0.46 - 0.12) + 39 * (0.10 - 0.04)), 0.02)
    ck("Zd019 nominal scenario", 21.0, _ad((0.46 - 0.12) + 31.3 * (0.06 - 0.05) * 0.6), 0.05)
    from scipy.optimize import brentq as _bq
    _hSI, _kbSI, _meSI = 1.054571817e-34, 1.380649e-23, 9.1093837e-31
    _TF2 = 2 * math.pi * _hSI ** 2 * (1e14 * 1e4) / (3 * _meSI * _kbSI)
    _T = _bq(lambda T: T * math.log(380 * _TF2 / T) - _TF2, 1.0, _TF2 * 0.9)
    ck("Zd019 BKT n2d=1e14 m=3", 231.0, _T, 0.02)
    ck("Zd019 Cu+diamond l12=0.5", 236.0, _tb(0.158, 223.0, 0.1, 1916.0, 0.5, 1916.0), 0.03)

    # lifetime
    kT = KB * 300
    ck("Ea 1yr", 1.261, kT * math.log(3.156e7 / (2 * math.pi * HBAR / 0.2)), 0.005)
    ck("Ea 1s", 0.815, kT * math.log(1 / (2 * math.pi * HBAR / 0.2)), 0.005)

    # ---- Zd020: 结构 → alpha^2F 的三段桥 ----
    _sp3 = _iu.spec_from_file_location("a2f", os.path.join(ROOT, "src", "alpha2f_bridge.py"))
    _a2f = _iu.module_from_spec(_sp3); _sp3.loader.exec_module(_a2f)
    _mo, _di, _df = _a2f.morse(), _a2f.diatomic(), _a2f.deform()
    ck("Phi C-C (N/m)", 627.0, _mo["Phi_N_per_m"], 0.01)
    ck("a C-C (1/A)", 2.33, _mo["a_per_A"], 0.01)
    ck("wlog diatomic", 0.841, _di["wlog"], 0.01)
    ck("D C-C (eV)", 5.78, _df["D_eV"], 0.01)
    ck("g C-C q=1/a (eV)", 0.063, _df["g_bond_eV"], 0.02)
    ck("g C-C q=1/R0 (eV)", 0.122, _df["g_at_R0_eV"], 0.02)
    ck("lam ratio (kF vs 1/a)", 9.57, (_df["g_kF_eV"] / _df["g_bond_eV"]) ** 2, 0.01)

    # ---- Zd016 §3.5: Migdal 适用区 ----
    _sp4 = _iu.spec_from_file_location("mz", os.path.join(ROOT, "src", "migdal_zone.py"))
    _mz = _iu.module_from_spec(_sp4); _sp4.loader.exec_module(_mz)
    _GT = {n: w / ef for n, w, ef, _ in _mz.TABLE}
    ck("Migdal diamond", 0.017, _GT["金刚石（C）"], 0.01)
    ck("Migdal MgB2", 0.053, _GT["MgB$_2$"], 0.01)
    ck("Migdal flat", 0.20, _GT["平带材料"], 0.01)
    ck("Migdal molecular", 0.27, _GT["分子晶体"], 0.01)

    # ---- Zd020 §5.5: 材料表 ----
    _mt = _json.load(open(os.path.join(ROOT, "data", "material_table.json")))
    _by = {r["material"]: r for r in _mt["bulk"]}
    ck("Phi diamond (omega_log)", 305.0, _by["金刚石 C"]["Phi_N_m"], 0.01)
    ck("n_eff diamond", 0.43, _by["金刚石 C"]["n_eff"], 0.02)
    ck("D_H W (eV)", 18.2, _by["W"]["D_H_eV"], 0.01)
    ck("D_eff W (eV)", 1.81, _by["W"]["D_eff_eV"], 0.02)
    ck("g_eff Nb (eV)", 0.031, _by["Nb"]["g_eff_eV"], 0.02)
    ck("g_I C-W (eV)", 0.061, _mt["interface"]["g_I_eV"], 0.02)
    ck("M_eff req (g=0.43)", 49.0, _mt["interface"]["M_eff_req_for_g43"], 0.02)

    # ---- Zd021: 形变势的通道分解 ----
    _dc = _json.load(open(os.path.join(ROOT, "data", "deform_channels.json")))
    ck("identity max err", 0.0, _dc["identity"]["max_err"], 1e-6)
    _ff = {f["n_e"]: f for f in _dc["filling_factor"]}
    ck("epsF/t at half filling", 0.0, _ff[1.0]["epsF_over_t"], 1e-3, 5e-3)
    _dm = {r["material"]: r for r in _dc["materials"]}
    ck("Cu Xi_d (eV)", 0.95, _dc["cu"]["D_d_implied_eV"], 0.02)
    ck("Cu D_hop (eV)", 5.58, _dc["cu"]["D_hop_eV"], 0.02)
    ck("C implied n_e (pure hopping)", 1.062, _dm["金刚石 C"]["n_e_if_hopping"], 0.01)

    # ---- Zd022: 两矩够不够（过滤版统计） ----
    _es = _json.load(open(os.path.join(ROOT, "data", "eliashberg_stats.json")))
    ck("Eliashberg n used", 672, _es["n_used"], 0.001)
    ck("Eliashberg n filtered", 293, _es["n_filtered"], 0.001)
    ck("shape effect median", 1.00, _es["shape_effect"]["median"], 0.01)
    ck("shape med lam 0.3-0.6", 0.996, _es["by_bin"]["0.3-0.6"]["median"], 0.01)
    ck("shape med lam 1.5-5", 1.51, _es["by_bin"]["1.5-5.0"]["median"], 0.02)
    ck("corr wpeak", 0.40, _es["corr_absdev_wpeak_over_wlog"], 0.02)
    _ec = _json.load(open(os.path.join(ROOT, "data", "eliashberg_clean.json")))
    ck("clean n", 499, _ec["n_clean"], 0.001)
    ck("clean filtered", 184, _ec["n_filtered"], 0.001)
    ck("clean shape median", 0.999, _ec["shape_effect"]["median"], 0.01)
    ck("clean lam 0.3-0.6", 1.005, _ec["by_bin"]["0.3-0.6"]["median"], 0.01)
    _src = _json.load(open("/Users/oygb/Downloads/Zero/_staging/材料库/jarvis_epc_data_figshare_1058.json"))
    _hi = [r for r in _src if r.get("lamb", 0) > 1.5]
    def _neg(r):
        x = r["a2F_original_x"]; y = r["a2F_original_y"]
        return any(xi < 0 and yi > 0 for xi, yi in zip(x, y))
    _hic = [r for r in _hi if r.get("stability") == "stable" and not _neg(r)]
    ck("n lam>1.5", 190, len(_hi), 0.01)
    ck("n lam>1.5 stable&clean", 3, len(_hic), 0.01)

    # ---- Zd023: slab 的强耦合端 ----
    _sl = _json.load(open(os.path.join(ROOT, "data", "slab_eliashberg.json")))
    _s3 = next(r for r in _sl["rows"] if r["kappa"] == 3.0 and r["D_I"] == 5.0)
    ck("slab AD k=3 DI=5", 213, _s3["Tc_AD"], 0.02)
    ck("slab full k=3 DI=5", 309, _s3["Tc_full"], 0.02)
    ck("slab ratio k=3 DI=5", 1.45, _s3["ratio"], 0.02)
    _s2 = next(r for r in _sl["rows"] if r["kappa"] == 2.0 and r["D_I"] == 10.0)
    ck("slab ratio k=2 DI=10", 2.71, _s2["ratio"], 0.02)

    # ---- Zd024: D_I 从结构来 ----
    _di = _json.load(open(os.path.join(ROOT, "data", "interface_di.json")))
    ck("n CC overlap", 3.73, _di["n_hom_check"]["C-C_overlap"], 0.02)
    ck("n WW overlap", 5.37, _di["n_hom_check"]["W-W_overlap"], 0.02)
    ck("n CW", 4.58, _di["n_overlap"], 0.02)
    ck("D_I geo", 8.9, _di["D_I_geo"], 0.03)
    _s9 = next(r for r in _di["slab"] if r["D_I"] == 9.0)
    ck("slab Tc_full DI=9", 615.0, _s9["Tc_full"], 0.03)

    # ---- Zd025: 五对界面与复扫 ----
    _is = _json.load(open(os.path.join(ROOT, "data", "interface_scan.json")))
    _pr = {p["pair"]: p for p in _is["pairs"]}
    ck("D_I C-W", 8.86, _pr["C–W"]["D_I"], 0.02)
    ck("D_I C-Hf", 12.14, _pr["C–Hf"]["D_I"], 0.02)
    _hf = next(r for r in _is["slab"] if r["pair"] == "C–Hf" and r["kappa"] == 2.0)
    ck("Tc_full C-Hf k=2", 850.0, _hf["Tc_full"], 0.03)
    _ws = _json.load(open(os.path.join(ROOT, "data", "window_scan.json")))
    ck("rescan JARVIS pass300", 0, _ws["JARVIS"]["n_pass_300"], 0.001)
    ck("rescan JARVIS max corr", 40.2, _ws["JARVIS"]["max_Tc_corr"], 0.02)
    ck("rescan MP pass300", 0, _ws["MP"]["n_pass_300"], 0.001)

    # ---- Zd026: bcc d 带 N_B ----
    _db = _json.load(open(os.path.join(ROOT, "data", "dband_dos.json")))
    _nb = {r["metal"]: r["N_B_per_eV"] for r in _db["rows"]}
    ck("N_B W", 1.536, _nb["W"], 0.02)
    ck("N_B Mo", 4.656, _nb["Mo"], 0.02)
    ck("N_B Hf", 0.665, _nb["Hf"], 0.02)
    ck("N_B D_I^2 Mo", 398.0, _nb["Mo"] * 9.24 ** 2, 0.02)

    # ---- Zd027: 界面配对的宇称 ----
    _pf = _json.load(open(os.path.join(ROOT, "data", "parity_form_factor.json")))
    ck("sp3 s weight", 0.25, _pf["sp3"]["w_s"], 0.01)
    ck("CM even weight", 0.25, _pf["CM_channel"]["sigma_pair_even_weight"], 0.01)
    ck("CM odd weight", 0.75, _pf["CM_channel"]["sigma_pair_odd_weight"], 0.01)

    # ---- Zd028: 超晶格几何 ----
    _sg = _json.load(open(os.path.join(ROOT, "data", "superlattice_geometry.json")))
    ck("xi 300K vF=1.33e6", 6.1, _sg["xi_nm"]["v_F=1.33e6"]["300"], 0.02)
    ck("xi 300K vF=4.6e5", 2.1, _sg["xi_nm"]["v_F=4.6e5"]["300"], 0.02)
    ck("f at d=20", 0.61, next(r["f"] for r in _sg["frac"] if r["d_nm"] == 20), 0.02)
    ck("d_max 300K", 12.2, _sg["window"]["d_max_nm"], 0.02)

    # ---- 家用条件：可电镀候选 ----
    _ep = _json.load(open(os.path.join(ROOT, "data", "electroplate_tc.json")))
    _cu = next(r for r in _ep["rows"] if r["metal"] == "Cu" and "金刚石" in r["mode"])
    ck("Cu D_I", 7.53, _cu["D_I"], 0.02)
    ck("Cu lam", 1.01, _cu["lam"], 0.02)
    ck("Cu Tc_AD diamond", 85, _cu["Tc_AD"], 0.03)
    _pt = next(r for r in _ep["rows"] if r["metal"] == "Pt" and "金刚石" in r["mode"])
    ck("Pt lam", 17.0, _pt["lam"], 0.03)

    # ---- 误差预算与界面模 ----
    _eb = _json.load(open(os.path.join(ROOT, "data", "tc_error_budget.json")))
    _mid = next(r for r in _eb["interface_mode"] if "中" in r["bond"])
    ck("omega_I Ag-C mid cm", 710, _mid["omega_cm"], 0.02)
    ck("omega_I Ag-C mid K", 1021, _mid["omega_K"], 0.02)
    ck("MC phys median", 118, _eb["bands"]["物理界面模（600–1200 K）"]["median"], 0.03)
    ck("sens D_I", 2.36, _eb["sensitivity"]["D_I"], 0.03)
    ck("sens mu*", -0.43, _eb["sensitivity"]["μ*"], 0.05)
    # ---- Ag/GO 全链条复算 ----
    _agrt = 0.1660 * (3.0 / (4 * 5.5)) * 8.6726 ** 2
    ck("Ag lambda(E_F=5.5)", 1.70, _agrt, 0.02)
    # ---- Zd029: 完整链 ----
    import math as _m
    def _ad(_l, _w, _mu=0.12):
        return _w / 1.2 * _m.exp(-1.04 * (1 + _l) / (_l - _mu * (1 + 0.62 * _l)))
    ck("Ag Tc_AD E_F=7", 97, _ad(0.1660 * (3 / 28) * 8.6726 ** 2, 1021.0), 0.02)
    ck("Ag Tc_AD E_F=5.5", 123, _ad(_agrt, 1021.0), 0.02)
    ck("Ag N_B floor", 0.107, 3.0 / 28.0, 0.01)

    # ---- Zd029 §4.5: FeSe/STO 反解 ----
    _fs = _json.load(open(os.path.join(ROOT, "data", "fese_benchmark.json")))
    _c65 = next(c for c in _fs["cases"] if c["case"] == "单层 FeSe/STO" and c["mu_star"] == 0.12)
    ck("FeSe/STO lambda 65K", 0.99, _c65["lambda_req"], 0.02)
    _c8 = next(c for c in _fs["cases"] if c["case"] == "体相 FeSe" and c["mu_star"] == 0.12)
    ck("FeSe bulk lambda 8K", 0.72, _c8["lambda_req"], 0.02)
    ck("Fe-O short-range lambda", 14.9, _fs["FeO_short_range"]["lambda"], 0.04)

    # ---- 双组份宇称平衡 ----
    _tp = _json.load(open(os.path.join(ROOT, "data", "two_component_parity.json")))
    _top = _tp["rows"][0]
    ck("two-comp top Tc", 658, _top["Tc"], 0.02)
    _s = next(r for r in _tp["rows"] if r["parity_ok"])
    ck("two-comp best parity Tc", 395, _s["Tc"], 0.02)
    ck("Debye base lam=0.5", 1.265, _es["debye_base"]["0.5"], 0.01)
    ck("Debye base lam=2.0", 1.294, _es["debye_base"]["2.0"], 0.01)
    ck("BCS self-test Tc", 2.70, _es["bcs_check"]["Tc_K"], 0.02)

    out = {"n_claims": len(claims), "n_pass": sum(c["pass"] for c in claims),
           "claims": claims}
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    for c in claims:
        print(f"{'PASS' if c['pass'] else 'FAIL'}  {c['claim']:28s} text={c['text']}  calc={c['calc']}")
    print(f"{out['n_pass']}/{out['n_claims']} PASS")
    print("saved", OUT)


if __name__ == "__main__":
    main()

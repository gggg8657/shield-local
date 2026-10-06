"""shield local 계산 엔진 — 결정론적 점선원 감마 차폐 계산 (stdlib 만, LLM 없음).

식 (점등방선원, 차폐는 선원–검출점 사이 평판, 공기 감쇠 무시):
  φ_i   = A·y_i / (4π r²)                                  비충돌 플루언스율 [cm⁻² s⁻¹]
  K̇_i   = φ_i · E_i · (μen/ρ)_air(E_i) · 1.602177e-13 J/MeV · 1000 g/kg · 3600 s/h   공기커마율 [Gy/h]
  X_i   = Σ_layers μ(E_i)·t                                 광학 두께 [mfp]  (μ = XCOM 간섭성 산란 제외 — ANS-6.4.3 규약)
  K̇     = Σ_i K̇_i · B(E_i, X_i) · exp(−X_i)
  Ḣ*(10) = Σ_i K̇_i,차폐 · (H*(10)/Ka)(E_i)                 ICRP 74 Table A.21 환산계수
데이터 출처는 data/*.json 의 _source 와 README 참고.
"""
import json
import math
import os

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MEV_J = 1.602176634e-13
CI = 3.7e10
LN2 = math.log(2)


def _load(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return json.load(f)


XCOM = _load("xcom.json")["materials"]
AIR = _load("air_muen.json")["table"]
NUC = _load("nuclides.json")
GP = _load("buildup_ans643.json")["materials"]

# ICRP Publication 74 (1996) Table A.21 / ICRU 57 — 단색 광자 H*(10)/Ka [Sv/Gy]
ICRP74_HK = [(0.010, 0.008), (0.015, 0.26), (0.020, 0.61), (0.030, 1.10), (0.040, 1.47), (0.050, 1.67), (0.060, 1.74),
             (0.080, 1.72), (0.100, 1.65), (0.150, 1.49), (0.200, 1.40), (0.300, 1.31), (0.400, 1.26), (0.500, 1.23),
             (0.600, 1.21), (0.800, 1.19), (1.0, 1.17), (1.5, 1.15), (2.0, 1.14), (3.0, 1.13), (4.0, 1.12), (5.0, 1.11),
             (6.0, 1.11), (8.0, 1.11), (10.0, 1.10)]

# 원자력안전법 시행령 별표 1 선량한도 (유효선량) — 개정 여부는 보건물리 담당자 확인
LIMITS = [
    {"id": "worker", "name": "방사선작업종사자", "annual_mSv": 50.0, "note": "연간 50 mSv 이하, 5년간 누적 100 mSv 이하(연평균 20 mSv)"},
    {"id": "visitor", "name": "수시출입자·운반종사자", "annual_mSv": 12.0, "note": "연간 12 mSv"},
    {"id": "public", "name": "일반인", "annual_mSv": 1.0, "note": "연간 1 mSv"},
]

ACT_UNITS = {"Bq": 1.0, "kBq": 1e3, "MBq": 1e6, "GBq": 1e9, "TBq": 1e12, "Ci": CI, "mCi": CI * 1e-3, "uCi": CI * 1e-6}
LEN_UNITS = {"mm": 0.1, "cm": 1.0, "m": 100.0, "in": 2.54}
RATE_UNITS = {"Sv/h": 1.0, "mSv/h": 1e-3, "uSv/h": 1e-6, "nSv/h": 1e-9, "Gy/h": 1.0, "mGy/h": 1e-3, "uGy/h": 1e-6, "nGy/h": 1e-9}
TIME_UNITS = {"s": 1 / 3600, "min": 1 / 60, "h": 1.0, "d": 24.0, "wk": 168.0, "y": 8766.0}
QUANT = {"kerma": ("공기커마율", "Gy/h"), "hstar": ("주위선량당량률 H*(10)", "Sv/h")}
RULES = {"max": "층별 재질 B(총 mfp) 중 최댓값 (보수적)", "last": "마지막 층 재질의 B(총 mfp)", "none": "축적인자 없음 (좁은 빔, 비보수적)"}


def norm_unit(u):
    return (u or "").strip().replace("µ", "u").replace("μ", "u").replace("㎝", "cm").replace("㎜", "mm")


def to_bq(value, unit):
    u = norm_unit(unit)
    if u not in ACT_UNITS:
        raise ValueError(f"방사능 단위를 모름: {unit} (Bq·kBq·MBq·GBq·TBq·Ci·mCi·uCi)")
    return float(value) * ACT_UNITS[u]


def to_cm(value, unit):
    u = norm_unit(unit)
    if u not in LEN_UNITS:
        raise ValueError(f"길이 단위를 모름: {unit} (mm·cm·m·in)")
    return float(value) * LEN_UNITS[u]


def rate_to_base(value, unit):
    """선량률 → Gy/h 또는 Sv/h (기본 단위)"""
    u = norm_unit(unit)
    if u not in RATE_UNITS:
        raise ValueError(f"선량률 단위를 모름: {unit}")
    return float(value) * RATE_UNITS[u]


def hours(value, unit):
    u = norm_unit(unit)
    if u not in TIME_UNITS:
        raise ValueError(f"시간 단위를 모름: {unit} (s·min·h·d·wk·y)")
    return float(value) * TIME_UNITS[u]


# ── 보간 ───────────────────────────────────────────────────────────────
def _loglog(x0, y0, x1, y1, x):
    return math.exp(math.log(y0) + (math.log(y1) - math.log(y0)) * (math.log(x) - math.log(x0)) / (math.log(x1) - math.log(x0)))


def _table(rows, e, col):
    """log-log 보간. 흡수끝(같은 E 두 줄)에서는 끝 아래 값(작은 μ, 보수적)을 쓴다."""
    if e < rows[0][0] * (1 - 1e-9) or e > rows[-1][0] * (1 + 1e-9):
        raise ValueError(f"에너지 {e * 1000:.1f} keV 가 표 범위({rows[0][0] * 1000:g} keV~{rows[-1][0]:g} MeV) 밖")
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a[0] <= e <= b[0] and b[0] > a[0]:
            return _loglog(a[0], a[col], b[0], b[col], e)
    return rows[-1][col]


def mu_rho(mat, e, coherent=False):
    """질량감쇠계수 μ/ρ [cm²/g] (NIST XCOM). coherent=False 가 차폐 계산용(ANS-6.4.3 규약)."""
    return _table(XCOM[mat]["table"], e, 1 if coherent else 2)


def muen_air(e):
    """공기 질량에너지흡수계수 μen/ρ [cm²/g] (NIST, Hubbell & Seltzer)"""
    return _table(AIR, e, 2)


def h_per_ka(e):
    """H*(10)/Ka [Sv/Gy] (ICRP 74). log-log 보간, 10 keV 미만은 0."""
    if e < ICRP74_HK[0][0]:
        return 0.0
    if e >= ICRP74_HK[-1][0]:
        return ICRP74_HK[-1][1]
    for (x0, y0), (x1, y1) in zip(ICRP74_HK, ICRP74_HK[1:]):
        if x0 <= e <= x1:
            return _loglog(x0, y0, x1, y1, e)


# ── G-P 축적인자 (ANSI/ANS-6.4.3, 노출) ─────────────────────────────────────
GP_XMAX = 40.0
_T2 = math.tanh(-2.0)


def gp_formula(b, c, a, xk, d, x):
    if x <= 0:
        return 1.0
    k = c * x ** a + d * (math.tanh(x / xk - 2.0) - _T2) / (1.0 - _T2)
    if abs(k - 1.0) < 1e-6:
        return 1.0 + (b - 1.0) * x
    return 1.0 + (b - 1.0) * (k ** x - 1.0) / (k - 1.0)


def _gp_at(t, j, x):
    return gp_formula(t["b"][j], t["c"][j], t["a"][j], t["Xk"][j], t["d"][j], x)


def buildup(gp_key, e, x):
    """B(E, x): 표 에너지 사이는 ln B 를 ln E 에 대해 선형 보간. x > 40 mfp 는 40 에서 자른다(경고)."""
    t = GP[gp_key]
    x = min(max(x, 0.0), GP_XMAX)
    es = t["E"]
    if e <= es[0]:
        return _gp_at(t, 0, x)
    if e >= es[-1]:
        return _gp_at(t, len(es) - 1, x)
    for j in range(len(es) - 1):
        if es[j] <= e <= es[j + 1]:
            b0, b1 = _gp_at(t, j, x), _gp_at(t, j + 1, x)
            return _loglog(es[j], b0, es[j + 1], b1, e)


# ── 선원 ───────────────────────────────────────────────────────────────
def nuclides():
    return {k: {"half_life_s": v["half_life_s"], "note": v["note"], "lines": len(v["photons"])} for k, v in NUC["gamma"].items()}


def photon_lines(nuclide=None, custom=None, cutoff_kev=20.0):
    """[(E MeV, 붕괴당 방출률, 종류, 출처)] — 차단에너지 미만 광자는 뺀다."""
    if custom:
        lines = [(float(e), float(y), "g", "사용자") for e, y in custom]
        for e, y, _, _ in lines:
            if not (0.01 <= e <= 10.0) or y <= 0:
                raise ValueError("사용자 정의 광자: 에너지 0.01~10 MeV, 방출률 > 0")
    else:
        if nuclide not in NUC["gamma"]:
            raise ValueError(f"지원하지 않는 핵종: {nuclide}")
        lines = [tuple(p) for p in NUC["gamma"][nuclide]["photons"]]
    return [p for p in lines if p[0] * 1000 >= cutoff_kev]


def half_life_s(nuclide):
    return NUC["gamma"].get(nuclide, {}).get("half_life_s")


# ── 차폐층 ─────────────────────────────────────────────────────────────
def norm_layers(layers):
    """[{material, thickness, unit, density?}] → [{mat, t_cm, rho}]"""
    out = []
    for l in layers or []:
        mat = l.get("material")
        if mat not in XCOM:
            raise ValueError(f"재질을 모름: {mat}")
        t = to_cm(l.get("thickness", 0), l.get("unit", "cm"))
        if t < 0:
            raise ValueError("두께는 0 이상")
        if t == 0:
            continue
        rho = float(l.get("density") or XCOM[mat]["rho"])
        out.append({"mat": mat, "t_cm": t, "rho": rho})
    return out


def _line_terms(e, layers, rule):
    """(X mfp, exp(−X), B, 적용 재질)"""
    x = sum(mu_rho(l["mat"], e) * l["rho"] * l["t_cm"] for l in layers)
    if not layers or rule == "none":
        return x, math.exp(-x), 1.0, None
    if rule == "last":
        g = XCOM[layers[-1]["mat"]]["buildup"]
        return x, math.exp(-x), buildup(g, e, x), g
    best = (1.0, None)
    for g in dict.fromkeys(XCOM[l["mat"]]["buildup"] for l in layers):
        b = buildup(g, e, x)
        if b > best[0]:
            best = (b, g)
    return x, math.exp(-x), best[0], best[1]


def point_rate(lines, a_bq, r_cm, layers, rule="max", detail=False):
    """{'kerma': Gy/h, 'hstar': Sv/h, 'lines': [...]} — layers 는 norm_layers 결과"""
    if r_cm <= 0:
        raise ValueError("거리는 0 보다 커야 함")
    geo = a_bq / (4 * math.pi * r_cm ** 2)
    tot_k = tot_h = 0.0
    rows = []
    for e, y, kind, origin in lines:
        k0 = geo * y * e * muen_air(e) * MEV_J * 1000.0 * 3600.0
        x, att, b, g = _line_terms(e, layers, rule)
        k = k0 * att * b
        h = k * h_per_ka(e)
        tot_k += k
        tot_h += h
        if detail:
            rows.append({"E": e, "y": y, "kind": kind, "origin": origin, "muen": muen_air(e), "k0": k0, "mfp": x,
                         "att": att, "B": b, "gp": g, "k": k, "hk": h_per_ka(e), "h": h, "kerma": k, "hstar": h})
    return {"kerma": tot_k, "hstar": tot_h, "lines": rows}


def gamma_constant(lines):
    """공기커마율 상수 Γ_δ [μGy·m²/(GBq·h)] — 비차폐, 1 m, 1 GBq"""
    return point_rate(lines, 1e9, 100.0, [])["kerma"] * 1e6


# ── 역산 ───────────────────────────────────────────────────────────────
def _bisect(f, target, lo, hi, it=200):
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) > target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-7 * max(hi, 1e-6):
            break
    return 0.5 * (lo + hi)


def solve_thickness(lines, a_bq, r_cm, fixed, mat, rho, target, quantity="hstar", rule="max", tmax_cm=2000.0):
    """fixed 층 뒤에 재질 mat 을 덧대 선량률이 target 이 되는 두께 [cm]. 도달 불가면 None."""
    if mat not in XCOM:
        raise ValueError(f"재질을 모름: {mat}")
    rho = float(rho or XCOM[mat]["rho"])

    def f(t):
        lay = fixed + ([{"mat": mat, "t_cm": t, "rho": rho}] if t > 0 else [])
        return point_rate(lines, a_bq, r_cm, lay, rule)[quantity]
    if f(0.0) <= target:
        return 0.0
    hi = 0.1
    while f(hi) > target:
        hi *= 2
        if hi > tmax_cm:
            return None
    return _bisect(f, target, hi / 2 if hi > 0.1 else 0.0, hi)


def hvl_tvl(lines, mat, rho=None, quantity="hstar", rule="max"):
    """첫 반가층·1/10가층 [cm] — 넓은 빔(축적인자 포함)과 좁은 빔(B=1). 거리 무관(같은 거리에서의 비)."""
    out = {}
    for name, ru in (("broad", rule if rule != "none" else "max"), ("narrow", "none")):
        base = point_rate(lines, 1.0, 100.0, [], ru)[quantity]
        for lab, frac in (("hvl", 0.5), ("tvl", 0.1)):
            out[f"{lab}_{name}"] = solve_thickness(lines, 1.0, 100.0, [], mat, rho, base * frac, quantity, ru)
    return out


# ── 시간·누적선량 ─────────────────────────────────────────────────────────
def cumulative(rate_per_h, t_h, half_s=None, decay=True):
    """누적선량 = ∫ Ḋ0·e^{−λt} dt  (decay=False 면 Ḋ0·T)"""
    if not decay or not half_s:
        return rate_per_h * t_h
    lam = LN2 / (half_s / 3600.0)
    return rate_per_h * (1 - math.exp(-lam * t_h)) / lam


# ── 베타 (간이 경고) ──────────────────────────────────────────────────────
BETA_MATS = {"PMMA (아크릴)": (1.19, 6.5), "알루미늄": (2.699, 13), "물·조직": (1.0, 7.4), "유리": (2.5, 11), "납": (11.35, 82)}


def beta_range_gcm2(emax):
    """Katz–Penfold (1952) 최대비정 [g/cm²]"""
    if emax <= 2.5:
        return 0.412 * emax ** (1.265 - 0.0954 * math.log(emax))
    return 0.530 * emax - 0.106


def beta_info(key):
    b = NUC["beta"].get(key)
    if not b:
        raise ValueError(f"베타 핵종을 모름: {key}")
    out = []
    for origin, q in b["Q_keV"]:
        emax = q / 1000.0
        r = beta_range_gcm2(emax)
        out.append({"origin": origin, "Emax_MeV": emax, "range_gcm2": r,
                    "thickness_mm": {m: 10 * r / rho for m, (rho, z) in BETA_MATS.items()},
                    "brems_frac": {m: 3.5e-4 * z * emax for m, (rho, z) in BETA_MATS.items()}})
    return {"key": key, "half_life_s": b["half_life_s"], "branches": out}


# ── 중성자 (비차폐 선량률만) ───────────────────────────────────────────────
# ISO 8529-1 기준 스펙트럼 평균 h*(10) [pSv·cm²] (ICRP 74 환산계수로 가중), 방출률 대표값
NEUTRON = {
    "Cf-252": {"h_pSv_cm2": 385.0, "yield": {"ug": 2.31e6, "Bq": 0.1165}, "note": "자발핵분열 분기 3.092 % × ν̄ 3.7676 = 붕괴당 0.1165 n → 2.31×10⁶ n/s/μg"},
    "Am-Be": {"h_pSv_cm2": 391.0, "yield": {"Ci": 2.2e6}, "note": "Am-241 1 Ci 당 약 2.2×10⁶ n/s — 선원 구조(Am:Be 비)에 따라 ±15 % 이상 다름. 인증서의 방출률을 쓰는 것이 좋다"},
}


def neutron_rate(source, amount, unit, r_cm):
    n = NEUTRON.get(source)
    if not n:
        raise ValueError(f"중성자 선원을 모름: {source}")
    u = norm_unit(unit)
    if u == "n/s":
        s = float(amount)
    elif u in ("ug", "mg", "g") and "ug" in n["yield"]:
        s = float(amount) * n["yield"]["ug"] * {"ug": 1, "mg": 1e3, "g": 1e6}[u]
    elif u in ACT_UNITS:
        bq = to_bq(amount, u)
        s = bq * n["yield"]["Bq"] if "Bq" in n["yield"] else bq / CI * n["yield"]["Ci"]
    else:
        raise ValueError("중성자 선원 양 단위: n/s, ug·mg(Cf-252), Bq·Ci 등")
    phi = s / (4 * math.pi * r_cm ** 2)
    return {"S_n_per_s": s, "phi": phi, "h_pSv_cm2": n["h_pSv_cm2"], "hstar": phi * n["h_pSv_cm2"] * 1e-12 * 3600.0, "note": n["note"]}


# ── 표시 ───────────────────────────────────────────────────────────────
def fmt_rate(v, quantity):
    """Gy/h·Sv/h 값을 읽기 쉬운 접두어로"""
    base = QUANT[quantity][1]
    s, rest = base[:2], base[2:]
    for f, p in ((1, ""), (1e-3, "m"), (1e-6, "μ"), (1e-9, "n")):
        if abs(v) >= f or f == 1e-9:
            return f"{v / f:.4g} {p}{s}{rest}"


def fmt_dose(v, quantity):
    return fmt_rate(v, quantity).replace("/h", "")


def fmt_len(cm):
    if cm is None:
        return "도달 불가"
    return f"{cm * 10:.3g} mm" if cm < 1 else f"{cm:.4g} cm"

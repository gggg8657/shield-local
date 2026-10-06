#!/usr/bin/env python3
"""shield local — 방사선 차폐·선량 간이 계산기. stdlib 만, 외부 전송 없음.

  python3 app.py                                   # http://localhost:8784
  LLM_API=ollama LLM_BASE_URL=http://localhost:11436 LLM_MODEL=gemma4:31b python3 app.py
  python3 app.py --cli '{"source":{"nuclide":"Cs-137","activity":{"value":10,"unit":"mCi"}},"distance":{"value":30,"unit":"cm"}}'

계산은 전부 shield.py (결정론). LLM 은 ① 자연어 질문 → 계산 입력 JSON ② 계산 결과 설명 문장에만 쓴다.
LLM 이 쓴 숫자는 계산하지 않으며, 설명문에 결과에 없는 숫자가 나오면 경고한다.
"""
import datetime
import json
import math
import os
import re
import secrets
import sys
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import shield as S
import validation

ROOT = os.path.dirname(os.path.abspath(__file__))
WS = os.environ.get("WORKSPACE") or os.path.join(ROOT, "_workspace")  # 포털이 AGENT_DATA/<도구> 로 모아 줌
LLM_API = os.environ.get("LLM_API", "ollama")
LLM_BASE = os.environ.get("LLM_BASE_URL", "http://localhost:8000/v1" if LLM_API == "openai" else "http://localhost:11434").rstrip("/")
MODEL = os.environ.get("LLM_MODEL", "qwen3:8b")
LLM_KEY = os.environ.get("LLM_API_KEY", "")
PORT = int(os.environ.get("PORT", "8784"))
LOCK = threading.Lock()
VERSION = "1.0"

DISCLAIMER = ("간이 계산 결과입니다. 점선원·평판 차폐·무한매질 축적인자 근사를 쓰며 선원 캡슐·산란 구조물·공기 감쇠는 반영하지 않습니다. "
              "인허가·정식 차폐 평가는 MCNP 등 정밀 해석과 보건물리 담당자 확인이 필요합니다.")


def read(p):
    with open(p, encoding="utf-8") as f:
        return f.read()


# ── 계산 ───────────────────────────────────────────────────────────────
def _num(d, key="value"):
    try:
        v = float((d or {}).get(key))
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _geom(r0, rmin, k=60, span=10.0):
    lo, hi = max(rmin, r0 / span), r0 * span
    return [lo * (hi / lo) ** (i / (k - 1)) for i in range(k)]


def calc(req):
    """POST /api/calc 본체 (mode=gamma). 결과에는 단계별 근거·그래프 데이터·경고가 들어 있다."""
    src = req.get("source") or {}
    q = req.get("quantity") if req.get("quantity") in S.QUANT else "hstar"
    rule = req.get("rule") if req.get("rule") in S.RULES else "max"
    cutoff = float(src.get("cutoff_kev") or 20.0)
    act = src.get("activity") or {}
    if _num(act) is None or _num(act) <= 0:
        raise ValueError("방사능 값을 입력하세요")
    a_bq = S.to_bq(_num(act), act.get("unit", "Bq"))
    custom = src.get("custom") or None
    nuc = None if custom else src.get("nuclide")
    lines = S.photon_lines(nuc, custom, cutoff)
    if not lines:
        raise ValueError("차단에너지 이상 광자선이 없습니다 (순수 베타 방출체는 '베타·중성자' 탭)")
    dist = req.get("distance") or {}
    if _num(dist) is None or _num(dist) <= 0:
        raise ValueError("거리를 입력하세요")
    r = S.to_cm(_num(dist), dist.get("unit", "cm"))
    layers = S.norm_layers(req.get("layers"))
    t_sum = sum(l["t_cm"] for l in layers)
    warns = []
    if t_sum >= r:
        raise ValueError(f"차폐 두께 합({S.fmt_len(t_sum)})이 선원–검출점 거리({S.fmt_len(r)}) 이상입니다")

    un = S.point_rate(lines, a_bq, r, [], rule, detail=True)
    sh = S.point_rate(lines, a_bq, r, layers, rule, detail=True)
    gam = S.gamma_constant(lines)
    half = None if custom else S.half_life_s(nuc)

    # 기여 큰 선 (차폐 후 기준)
    order = sorted(range(len(sh["lines"])), key=lambda i: -sh["lines"][i][q])
    tot = sh[q] or 1e-300
    line_rows = []
    for i in order[:40]:
        a, b = un["lines"][i], sh["lines"][i]
        line_rows.append({"E_keV": a["E"] * 1000, "y": a["y"], "kind": a["kind"], "origin": a["origin"], "k0": a["k0"], "h0": a["h"],
                          "mfp": b["mfp"], "B": b["B"], "gp": b["gp"], "k": b["k"], "h": b["h"], "share": b[q] / tot})

    # 경고
    if any(b["mfp"] > S.GP_XMAX for b in sh["lines"] if b[q] > 1e-4 * tot):
        warns.append(f"일부 광자선의 광학 두께가 {S.GP_XMAX:g} mfp 를 넘어 축적인자를 {S.GP_XMAX:g} mfp 값으로 잘랐습니다(그 선의 기여는 극히 작음).")
    low = sum(b[q] for b in sh["lines"] if b["E"] < 0.03)
    if low > 0.2 * tot:
        warns.append(f"30 keV 미만 저에너지 광자가 선량률의 {100 * low / tot:.0f}% 를 차지합니다. 선원 캡슐·용기 벽의 자체 감쇠가 크게 작용하므로 실제 값은 훨씬 낮을 수 있습니다(차단에너지 설정 참고).")
    subs = {l["mat"] for l in layers if l["mat"] in ("polyethylene", "ss304")}
    if subs:
        warns.append("축적인자 대용: " + ", ".join({"polyethylene": "폴리에틸렌 → 물의 G-P 계수", "ss304": "SS304 → 철의 G-P 계수"}[m] for m in sorted(subs)) + " (ANS-6.4.3 에 해당 재질 없음).")
    if len({l["mat"] for l in layers}) > 1 and rule != "none":
        warns.append(f"다층 차폐: 축적인자는 '{S.RULES[rule]}' 근사입니다. 저Z 재질이 마지막 층이면 B 가 크게 잡힙니다.")
    if rule == "none" and layers:
        warns.append("축적인자를 끈 좁은 빔 계산은 산란선을 빼므로 선량률을 과소평가합니다.")
    if q == "hstar" and any(l["mat"] in ("concrete", "water", "polyethylene", "aluminium") for l in layers):
        warns.append("H*(10) 은 원래 광자 에너지의 환산계수(ICRP 74)를 산란선까지 포함한 공기커마에 곱한 근사입니다. 두꺼운 저Z 차폐에서는 산란선 저에너지 성분 때문에 다소 과소평가될 수 있습니다.")
    if any(b["E"] < 0.03 and b["gp"] == "lead" for b in sh["lines"]) or any(0.08 < b["E"] < 0.1 for b in sh["lines"] if b["gp"] == "lead" and b[q] > 0.05 * tot):
        warns.append("납 K 흡수끝(88 keV) 부근·30 keV 미만 광자의 납 축적인자는 표 경계값 보간입니다.")

    res = {"version": VERSION, "quantity": q, "quantity_name": S.QUANT[q][0], "unit": S.QUANT[q][1], "rule": rule, "rule_name": S.RULES[rule],
           "nuclide": nuc or "사용자 정의", "half_life_s": half, "note": (S.NUC["gamma"].get(nuc) or {}).get("note", ""),
           "activity_bq": a_bq, "activity_in": act, "r_cm": r, "distance_in": dist, "cutoff_kev": cutoff,
           "layers": [dict(l, name=S.XCOM[l["mat"]]["name"]) for l in layers], "n_lines": len(lines),
           "gamma_const": gam, "unshielded": {"kerma": un["kerma"], "hstar": un["hstar"]}, "shielded": {"kerma": sh["kerma"], "hstar": sh["hstar"]},
           "transmission": sh[q] / un[q] if un[q] else None, "lines": line_rows, "warnings": warns, "disclaimer": DISCLAIMER}

    # 단계별 근거 (기여 상위 3개 선)
    res["steps"] = steps(res, lines, a_bq, r, layers, rule, q, [i for i in order[:3] if i == order[0] or sh["lines"][i][q] >= 1e-3 * tot], un, sh)

    # 역산: 목표 선량률
    tgt = req.get("target") or {}
    smat = req.get("solve_material") or (layers[-1]["mat"] if layers else "lead")
    srho = float(req.get("solve_density") or S.XCOM[smat]["rho"])
    if _num(tgt) is not None and _num(tgt) > 0:
        tv = S.rate_to_base(_num(tgt), tgt.get("unit", "uSv/h"))
        t_need = S.solve_thickness(lines, a_bq, r, layers, smat, srho, tv, q, rule)
        need_r = r * math.sqrt(sh[q] / tv) if sh[q] > tv else None
        res["target"] = {"value": tv, "in": tgt, "material": smat, "material_name": S.XCOM[smat]["name"], "rho": srho, "t_cm": t_need,
                         "r_cm_for_target": need_r, "already": sh[q] <= tv}
        if t_need is not None and t_sum + t_need >= r:
            warns.append("필요 두께가 선원–검출점 거리보다 큽니다. 거리를 늘리거나 선원 쪽 차폐를 검토하세요.")
    hv = S.hvl_tvl(lines, smat, srho, q, rule)
    deep = [S.solve_thickness(lines, 1.0, 100.0, [], smat, srho, S.point_rate(lines, 1.0, 100.0, [], rule)[q] * f, q, rule) for f in (1e-2, 1e-3)]
    hv["tvl_eq"] = (deep[1] - deep[0]) if None not in deep else None
    res["hvl"] = dict(hv, material=smat, material_name=S.XCOM[smat]["name"], rho=srho)

    # 시간·누적선량 (선량한도 비교는 H*(10))
    tm = req.get("time") or {}
    if _num(tm) is not None and _num(tm) > 0:
        th = S.hours(_num(tm), tm.get("unit", "h"))
        decay = bool(req.get("decay", True))
        dose = {k: S.cumulative(sh[k], th, half, decay) for k in ("kerma", "hstar")}
        res["time"] = {"hours": th, "in": tm, "decay": decay and bool(half), "dose": dose,
                       "limits": [dict(l, frac=dose["hstar"] * 1e3 / l["annual_mSv"]) for l in S.LIMITS]}
        if decay and half and th > 0.05 * half / 3600:
            res["time"]["decay_note"] = f"반감기 {fmt_time(half)} — 작업시간 동안 붕괴를 반영함 (붕괴 무시 시 {S.fmt_dose(sh['hstar'] * th, 'hstar')})"

    # 그래프 데이터
    tmax = max([x for x in (res.get("target", {}).get("t_cm"), (hv.get("tvl_broad") or 1) * 3) if x] + [0.5])
    tmax = min(tmax * 1.3, max(r - t_sum, 0.1) * 0.999)
    res["curve_t"] = {"material": smat, "material_name": S.XCOM[smat]["name"], "points": [
        [t, S.point_rate(lines, a_bq, r, layers + [{"mat": smat, "t_cm": t, "rho": srho}] if t > 0 else layers, rule)[q]]
        for t in [tmax * i / 60 for i in range(61)]]}
    rs = _geom(r, t_sum * 1.001 + 0.01)
    res["curve_r"] = {"points": [[x, sh[q] * (r / x) ** 2] for x in rs]}
    return res


def steps(res, lines, a_bq, r, layers, rule, q, idx, un, sh):
    act = res["activity_in"]
    out = [{"t": "방사능 환산", "f": "A = 값 × 단위",
            "v": f"A = {act.get('value')} {act.get('unit')} = {a_bq:.4g} Bq" + ("  (1 Ci = 3.7×10¹⁰ Bq)" if "Ci" in str(act.get("unit")) else "")},
           {"t": "거리", "f": "r", "v": f"r = {r:.4g} cm"}]
    for n, i in enumerate(idx, 1):
        a, b = un["lines"][i], sh["lines"][i]
        e = a["E"]
        phi = a_bq * a["y"] / (4 * math.pi * r * r)
        out.append({"t": f"[선 {n}] {e * 1000:.2f} keV ({a['origin']}, 방출률 {a['y'] * 100:.4g} %)",
                    "f": "φ = A·y / (4π r²)", "v": f"φ = {a_bq:.4g} × {a['y']:.5g} / (4π × {r:.4g}²) = {phi:.4g} cm⁻² s⁻¹"})
        out.append({"t": "", "f": "K̇₀ = φ · E · (μen/ρ)공기 · 1.602×10⁻¹³ J/MeV · 10³ g/kg · 3600 s/h",
                    "v": f"(μen/ρ)공기({e * 1000:.1f} keV) = {a['muen']:.5g} cm²/g (NIST) → K̇₀ = {S.fmt_rate(a['k0'], 'kerma')}"})
        if layers:
            parts = [f"{S.mu_rho(l['mat'], e):.5g}×{l['rho']:g}×{l['t_cm']:.4g}" for l in layers]
            out.append({"t": "", "f": "X = Σ (μ/ρ)·ρ·t   (μ/ρ: NIST XCOM, 간섭성 산란 제외)",
                        "v": f"X = {' + '.join(parts)} = {b['mfp']:.4g} mfp → e^(−X) = {b['att']:.4g}"})
            out.append({"t": "", "f": f"B = G-P(ANS-6.4.3, {b['gp'] or '-'}; E, X) — {S.RULES[rule]}",
                        "v": f"B = {b['B']:.4g} → K̇ = K̇₀ · B · e^(−X) = {S.fmt_rate(b['k'], 'kerma')}"})
        out.append({"t": "", "f": "Ḣ*(10) = K̇ · (H*(10)/Ka)(E)  (ICRP 74)",
                    "v": f"H*(10)/Ka = {b['hk']:.4g} Sv/Gy → Ḣ*(10) = {S.fmt_rate(b['h'], 'hstar')}"})
    rest = len(lines) - len(idx)
    out.append({"t": "합계", "f": f"모든 광자선 {len(lines)}개 합" + (f" (위 {len(idx)}개 외 {rest}개 포함)" if rest > 0 else ""),
                "v": f"비차폐 {S.fmt_rate(un[q], q)} → 차폐 후 {S.fmt_rate(sh[q], q)}"
                     + (f" (투과율 {sh[q] / un[q]:.4g})" if layers and un[q] else "")})
    return out


def fmt_time(s):
    for f, u in ((365.25 * 86400, "년"), (86400, "일"), (3600, "시간"), (60, "분")):
        if s >= f:
            return f"{s / f:.4g} {u}"
    return f"{s:.3g} 초"


def calc_other(req):
    if req.get("mode") == "beta":
        return {"mode": "beta", **S.beta_info(req.get("nuclide")), "disclaimer": DISCLAIMER}
    d = req.get("distance") or {}
    r = S.to_cm(_num(d) or 100, d.get("unit", "cm"))
    a = req.get("amount") or {}
    out = S.neutron_rate(req.get("source"), _num(a), a.get("unit"), r)
    tm = req.get("time") or {}
    if _num(tm):
        out["dose"] = out["hstar"] * S.hours(_num(tm), tm.get("unit", "h"))
        out["limits"] = [dict(l, frac=out["dose"] * 1e3 / l["annual_mSv"]) for l in S.LIMITS]
    return {"mode": "neutron", "source": req.get("source"), "r_cm": r, **out, "disclaimer": DISCLAIMER,
            "shielding": "중성자 차폐 계산은 지원하지 않습니다(제거단면 근사는 오차가 커서 제외). 정밀 해석이 필요합니다."}


# ── 보고서 (Markdown) ────────────────────────────────────────────────────
def report_md(res, title=""):
    q = res["quantity"]
    L = [f"# 방사선 차폐 간이 계산 보고서{(' — ' + title) if title else ''}", "",
         f"- 작성: {datetime.datetime.now():%Y-%m-%d %H:%M} · shield local v{res.get('version', VERSION)}",
         f"- **주의**: {res['disclaimer']}", "", "## 입력", "",
         f"| 항목 | 값 |", "|---|---|",
         f"| 선원 | {res['nuclide']}{(' (' + res['note'] + ')') if res.get('note') else ''} |",
         f"| 방사능 | {res['activity_in'].get('value')} {res['activity_in'].get('unit')} = {res['activity_bq']:.4g} Bq |",
         f"| 거리 | {res['distance_in'].get('value')} {res['distance_in'].get('unit')} ({res['r_cm']:.4g} cm) |",
         "| 차폐 | " + (" → ".join(f"{l['name']} {S.fmt_len(l['t_cm'])} (ρ={l['rho']:g})" for l in res["layers"]) or "없음") + " |",
         f"| 선량 양 | {res['quantity_name']} |", f"| 축적인자 | {res['rule_name']} |",
         f"| 광자 차단에너지 | {res['cutoff_kev']:g} keV (광자선 {res['n_lines']}개) |", "",
         "## 결과", "",
         f"- 비차폐: 공기커마율 **{S.fmt_rate(res['unshielded']['kerma'], 'kerma')}**, H*(10) **{S.fmt_rate(res['unshielded']['hstar'], 'hstar')}**",
         f"- 차폐 후: 공기커마율 **{S.fmt_rate(res['shielded']['kerma'], 'kerma')}**, H*(10) **{S.fmt_rate(res['shielded']['hstar'], 'hstar')}**",
         f"- 공기커마율 상수 Γ(δ={res['cutoff_kev']:g} keV) = {res['gamma_const']:.4g} μGy·m²/(GBq·h)"]
    t = res.get("target")
    if t:
        L.append(f"- 목표 {S.fmt_rate(t['value'], q)} → {t['material_name']}(ρ={t['rho']:g}) 추가 두께 **{S.fmt_len(t['t_cm'])}**"
                 + (f"; 차폐 그대로 거리만 늘리면 **{S.fmt_len(t['r_cm_for_target'])}**" if t.get("r_cm_for_target") else " (이미 목표 이하)"))
    h = res.get("hvl")
    if h:
        L.append(f"- {h['material_name']} 반가층/1/10가층 (넓은 빔·첫 층): {S.fmt_len(h['hvl_broad'])} / {S.fmt_len(h['tvl_broad'])}; "
                 f"평형 TVL {S.fmt_len(h.get('tvl_eq'))}; 좁은 빔: {S.fmt_len(h['hvl_narrow'])} / {S.fmt_len(h['tvl_narrow'])}")
    tm = res.get("time")
    if tm:
        L.append(f"- 작업 {tm['in'].get('value')} {tm['in'].get('unit')}: 누적 H*(10) **{S.fmt_dose(tm['dose']['hstar'], 'hstar')}**"
                 + (" (붕괴 반영)" if tm.get("decay") else ""))
        for l in tm["limits"]:
            L.append(f"  - {l['name']} 연간 {l['annual_mSv']:g} mSv 의 {100 * l['frac']:.3g}%")
    L += ["", "## 계산 근거", ""]
    for s in res["steps"]:
        L.append(f"- {('**' + s['t'] + '** ') if s['t'] else ''}`{s['f']}`  \n  {s['v']}")
    L += ["", "## 광자선별 기여 (상위)", "", "| E (keV) | 방출률 | mfp | B | 차폐 후 | 비율 |", "|---:|---:|---:|---:|---:|---:|"]
    for x in res["lines"][:15]:
        L.append(f"| {x['E_keV']:.2f} | {x['y']:.4g} | {x['mfp']:.3g} | {x['B']:.3g} | {S.fmt_rate(x['h' if q == 'hstar' else 'k'], q)} | {100 * x['share']:.1f}% |")
    if res["warnings"]:
        L += ["", "## 근사·경고", ""] + [f"- {w}" for w in res["warnings"]]
    L += ["", "## 데이터 출처", "", *[f"- {s}" for s in SOURCES], ""]
    return "\n".join(L)


SOURCES = [
    "감쇠계수 μ/ρ: NIST XCOM v3.1 (Berger et al., NIST SRD 8) — 원소별 값, 혼합물은 NIST 조성(NISTIR 5632 Table 2) 질량분율 가중합",
    "공기 μen/ρ: Hubbell & Seltzer, NISTIR 5632 (NIST X-Ray Mass Attenuation Coefficients, Air Dry)",
    "축적인자: ANSI/ANS-6.4.3-1991 G-P 노출 축적인자 계수 (NUREG/CR-5740, Trubey 1991, Table 5.1)",
    "H*(10)/Ka: ICRP Publication 74 (1996) Table A.21",
    "붕괴 데이터: DDEP 권고값 (LNHB/CEA LARA)",
    "선량한도: 원자력안전법 시행령 별표 1",
]


# ── 이력 ───────────────────────────────────────────────────────────────
def hist_path(hid):
    if not re.fullmatch(r"\d{8}-\d{6}-[0-9a-f]{4}", hid or ""):
        raise ValueError("잘못된 이력 id")
    return os.path.join(WS, "history", hid + ".json")


def hist_save(kind, req, res, title):
    os.makedirs(os.path.join(WS, "history"), exist_ok=True)
    now = datetime.datetime.now()
    rec = {"id": f"{now:%Y%m%d-%H%M%S}-{secrets.token_hex(2)}", "ts": now.isoformat(timespec="seconds"), "kind": kind,
           "title": title, "request": req, "result": res}
    p = hist_path(rec["id"])
    with LOCK:
        with open(p + ".tmp", "w", encoding="utf-8") as f:
            json.dump(rec, f, ensure_ascii=False)
        os.replace(p + ".tmp", p)
    return rec["id"]


def hist_list(limit=300):
    d = os.path.join(WS, "history")
    if not os.path.isdir(d):
        return []
    out = []
    for fn in sorted(os.listdir(d), reverse=True)[:limit]:
        if fn.endswith(".json"):
            try:
                r = json.loads(read(os.path.join(d, fn)))
                out.append({"id": r["id"], "ts": r["ts"], "kind": r["kind"], "title": r["title"]})
            except Exception:
                pass
    return out


def title_of(res):
    lay = " + ".join(f"{l['name'].split(' ')[0]} {S.fmt_len(l['t_cm'])}" for l in res["layers"]) or "비차폐"
    a = res["activity_in"]
    return f"{res['nuclide']} {a.get('value')} {a.get('unit')} · {S.fmt_len(res['r_cm'])} · {lay} → {S.fmt_rate(res['shielded'][res['quantity']], res['quantity'])}"


# ── LLM (보조: 질문 → 입력 JSON, 결과 → 설명) ─────────────────────────────────
def llm(system, user, temperature=0.0, on_token=lambda t: None, json_mode=False):
    msgs = [{"role": "system", "content": system}, {"role": "user", "content": user}]
    try:
        if LLM_API == "openai":
            body = {"model": MODEL, "stream": True, "temperature": temperature, "messages": msgs}
            hdr = {"Content-Type": "application/json", **({"Authorization": f"Bearer {LLM_KEY}"} if LLM_KEY else {})}
            buf = []
            with urllib.request.urlopen(urllib.request.Request(LLM_BASE + "/chat/completions", json.dumps(body).encode(), hdr), timeout=600) as r:
                for line in r:
                    line = line.decode().strip()
                    if line.startswith("data:") and line != "data: [DONE]":
                        tok = (json.loads(line[5:])["choices"][0].get("delta") or {}).get("content") or ""
                        if tok:
                            buf.append(tok)
                            on_token(tok)
            return "".join(buf)
        body = {"model": MODEL, "stream": True, "think": False, "messages": msgs, "options": {"temperature": temperature, "num_ctx": 8192}}
        if json_mode:
            body["format"] = "json"
        for attempt in (0, 1):
            try:
                buf = []
                req = urllib.request.Request(LLM_BASE + "/api/chat", json.dumps(body).encode(), {"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=600) as r:
                    for line in r:
                        if not line.strip():
                            continue
                        j = json.loads(line)
                        if "error" in j:
                            raise RuntimeError(j["error"])
                        tok = j.get("message", {}).get("content", "")
                        if tok:
                            buf.append(tok)
                            on_token(tok)
                        if j.get("done"):
                            break
                return "".join(buf)
            except urllib.error.HTTPError as e:
                msg = e.read().decode(errors="replace")
                if attempt == 0 and "think" in msg:
                    body.pop("think")
                    continue
                raise RuntimeError(f"LLM HTTP {e.code}: {msg[:300]}")
    except urllib.error.URLError as e:
        raise RuntimeError(f"LLM 서버 연결 실패 ({LLM_BASE}): {e.reason}")


def llm_ok():
    try:
        url = LLM_BASE + ("/models" if LLM_API == "openai" else "/api/tags")
        with urllib.request.urlopen(url, timeout=2):
            return True
    except Exception:
        return False


PARSE_SYS = """너는 방사선 차폐 계산기의 입력 변환기다. 사용자의 한국어 질문을 계산 입력 JSON 으로 바꾼다. 계산은 절대 하지 않는다.
JSON 하나만 출력한다. 키:
{"nuclide": 핵종 키(아래 목록 중 하나) 또는 null,
 "activity": {"value": 숫자, "unit": "Bq|kBq|MBq|GBq|TBq|Ci|mCi|uCi"} 또는 null,
 "distance": {"value": 숫자, "unit": "mm|cm|m"} 또는 null,
 "layers": [{"material": 재질 키, "thickness": 숫자, "unit": "mm|cm|m"}]  (질문에 이미 있는 차폐. 없으면 []),
 "target": {"value": 숫자, "unit": "uSv/h|mSv/h|nSv/h|Sv/h|uGy/h|mGy/h|Gy/h"} 또는 null,
 "solve_material": 두께를 구할 재질 키 또는 null,
 "quantity": "hstar"(선량당량·Sv, 기본) 또는 "kerma"(공기커마·Gy),
 "time": {"value": 숫자, "unit": "min|h|d|wk|y"} 또는 null,
 "ask": "thickness"(필요 두께) | "rate"(선량률) | "distance"(필요 거리) | "dose"(누적선량),
 "missing": [빠진 필수 정보 설명 문자열]}
규칙:
- 질문에 적힌 숫자와 단위를 그대로 옮긴다. 단위를 바꾸거나 숫자를 계산·추정하지 않는다. µ·μ 는 u 로 쓴다.
- "납 몇 mm/두께는?" 처럼 두께를 묻고 목표 선량률이 없으면 target 은 null 로 두고 missing 에 "목표 선량률" 을 넣는다.
- 핵종 이름은 목록의 키로 정규화한다(세슘-137, Cs137, 137Cs → "Cs-137"; 코발트 → "Co-60" 은 질량수가 없으면 missing 에 적는다).
- 재질: 납→lead, 텅스텐→tungsten, 철·강철→iron, 스테인리스→ss304, 구리→copper, 알루미늄→aluminium, 콘크리트→concrete, 물→water, 폴리에틸렌·PE→polyethylene.
- 필수: nuclide, activity, distance. 없으면 null 로 두고 missing 에 적는다."""


def parse_question(text, on_token=lambda t: None):
    nucs = ", ".join(S.NUC["gamma"])
    raw = llm(PARSE_SYS + "\n핵종 키 목록: " + nucs + "\n재질 키 목록: " + ", ".join(S.XCOM), "질문: " + text.strip(),
              on_token=on_token, json_mode=True)
    m = re.search(r"\{.*\}", re.sub(r"<think>.*?</think>", "", raw, flags=re.S), re.S)
    if not m:
        raise RuntimeError("LLM 이 JSON 을 돌려주지 않았습니다: " + raw[:200])
    return validate_parsed(json.loads(m.group(0)), text), raw


def _canon_nuc(n):
    if not n:
        return None
    s = str(n).strip().replace(" ", "")
    for k in S.NUC["gamma"]:
        el, a = k.split("-", 1)
        if s.lower() in (k.lower(), (el + a).lower(), (a + el).lower()):
            return k
    return None


def validate_parsed(p, text=""):
    """LLM 출력 → 계산 요청. 단위·키를 코드가 다시 검사하고, 숫자가 질문에 실제로 있는지 확인한다."""
    issues = list(p.get("missing") or [])
    nums_in_q = {float(x) for x in re.findall(r"\d+(?:\.\d+)?", text.replace(",", ""))}
    req = {"source": {}, "layers": [], "quantity": p.get("quantity") if p.get("quantity") in S.QUANT else "hstar"}
    nuc = _canon_nuc(p.get("nuclide"))
    if p.get("nuclide") and not nuc:
        issues.append(f"지원하지 않는 핵종: {p.get('nuclide')}")
    req["source"]["nuclide"] = nuc

    def chk(obj, units, label):
        if not isinstance(obj, dict) or obj.get("value") in (None, ""):
            return None
        u = S.norm_unit(obj.get("unit"))
        try:
            v = float(obj["value"])
        except (TypeError, ValueError):
            issues.append(f"{label} 값이 숫자가 아님")
            return None
        if u not in units:
            issues.append(f"{label} 단위를 모름: {obj.get('unit')}")
            return None
        if text and v not in nums_in_q:
            issues.append(f"{label} 값 {v:g} 이(가) 질문에 없음 — LLM 이 만든 숫자일 수 있어 확인 필요")
        return {"value": v, "unit": u}
    req["source"]["activity"] = chk(p.get("activity"), S.ACT_UNITS, "방사능")
    req["distance"] = chk(p.get("distance"), S.LEN_UNITS, "거리")
    for l in p.get("layers") or []:
        if l.get("material") in S.XCOM:
            t = chk({"value": l.get("thickness"), "unit": l.get("unit")}, S.LEN_UNITS, "차폐 두께")
            if t:
                req["layers"].append({"material": l["material"], "thickness": t["value"], "unit": t["unit"]})
        else:
            issues.append(f"재질을 모름: {l.get('material')}")
    tg = chk(p.get("target"), S.RATE_UNITS, "목표 선량률")
    if tg:
        req["target"] = tg
        if tg["unit"].endswith("Gy/h"):
            req["quantity"] = "kerma"
        elif tg["unit"].endswith("Sv/h"):
            req["quantity"] = "hstar"
    if p.get("solve_material") in S.XCOM:
        req["solve_material"] = p["solve_material"]
    elif p.get("solve_material"):
        issues.append(f"재질을 모름: {p.get('solve_material')}")
    tm = chk(p.get("time"), S.TIME_UNITS, "시간")
    if tm:
        req["time"] = tm
    for k, lab in (("nuclide", "핵종"), ("activity", "방사능")):
        if not req["source"].get(k):
            issues.append(f"{lab} 정보 필요")
    if not req.get("distance"):
        issues.append("거리 정보 필요")
    req["ask"] = p.get("ask") or "rate"
    return {"request": req, "issues": list(dict.fromkeys(issues))}


EXPLAIN_SYS = """너는 보건물리 계산 결과를 연구자에게 설명하는 조수다. 아래 [결과 요약]에 있는 숫자만 인용해 한국어로 4~7문장으로 설명한다.
- 숫자를 새로 계산하거나 반올림을 바꾸거나 단위를 바꾸지 않는다. 요약에 없는 숫자는 쓰지 않는다.
- 질문에 대한 답(필요 두께·선량률·거리·누적선량)을 첫 문장에 둔다. 그다음 근거(감쇠·축적인자)와 주요 근사·주의점을 짧게.
- 계산은 광자선별 지수감쇠 × 축적인자로 했다. 반가층·1/10가층은 참고값일 뿐이므로 "적용했다"고 쓰지 않는다.
- 마지막 문장은 항상 "간이 계산이므로 정식 평가는 정밀 해석과 보건물리 담당자 확인이 필요합니다."
- 마크다운 제목·표·이모지 금지."""


def summary_for_llm(res, ask):
    q = res["quantity"]
    s = [f"선원: {res['nuclide']} {res['activity_in'].get('value')} {res['activity_in'].get('unit')}",
         f"거리: {S.fmt_len(res['r_cm'])}",
         f"현재 차폐: {', '.join(l['name'] + ' ' + S.fmt_len(l['t_cm']) for l in res['layers']) or '없음'}",
         f"선량 양: {res['quantity_name']}",
         f"비차폐 선량률: {S.fmt_rate(res['unshielded'][q], q)}",
         f"현재 차폐 후 선량률: {S.fmt_rate(res['shielded'][q], q)}"]
    t = res.get("target")
    if t:
        s.append(f"목표 선량률: {S.fmt_rate(t['value'], q)}")
        s.append(f"목표를 위한 {t['material_name']} 추가 두께: {S.fmt_len(t['t_cm'])}")
        if t.get("r_cm_for_target"):
            s.append(f"차폐 없이 거리로만 맞추면 필요한 거리: {S.fmt_len(t['r_cm_for_target'])}")
    h = res.get("hvl")
    if h:
        s.append(f"참고값(계산에 쓰이지 않음) — {h['material_name']} 첫 반가층(넓은 빔): {S.fmt_len(h['hvl_broad'])}, 첫 1/10가층: {S.fmt_len(h['tvl_broad'])}")
    tm = res.get("time")
    if tm:
        s.append(f"작업시간 {tm['in'].get('value')} {tm['in'].get('unit')} 누적선량: {S.fmt_dose(tm['dose']['hstar'], 'hstar')}")
    if res["layers"]:
        s.append(f"축적인자 근사: {res['rule_name']}")
    else:
        s.append("차폐 없음: 거리 역제곱만 적용")
    s += [f"주의: {w}" for w in res["warnings"][:3]]
    return "\n".join(s) + f"\n질문 유형: {ask}"


def numbers_of(text):
    return {re.sub(r"[^\d.]", "", m).strip(".") for m in re.findall(r"\d[\d,]*(?:\.\d+)?", text or "")} - {""}


def ask(req, emit):
    """자연어 질문: 변환(LLM) → 검증(코드) → 계산(결정론) → 설명(LLM) + 숫자 검사"""
    text = (req.get("question") or "").strip()
    if not text:
        raise ValueError("질문을 입력하세요")
    emit({"stage": "parse"})
    parsed, raw = parse_question(text, lambda t: emit({"token": t, "ch": "parse"}))
    emit({"parsed": parsed, "raw": raw})
    blocking = [i for i in parsed["issues"] if "필요" in i or "모름" in i or "지원하지" in i]  # 질문에 없는 숫자(LLM 날조 의심)도 멈춘다
    creq = parsed["request"]
    if blocking or (creq.get("ask") == "thickness" and not creq.get("target")):
        return {"parsed": parsed, "result": None, "explain": None,
                "message": "계산에 필요한 정보가 부족하거나 확인이 필요합니다. '계산 탭에 채우기'로 옮겨 값을 확인·보완한 뒤 계산하세요."}
    emit({"stage": "calc"})
    res = calc(creq)
    hid = hist_save("ask", dict(creq, question=text), res, "Q: " + text[:60])
    emit({"result": res, "id": hid})
    emit({"stage": "explain"})
    summ = summary_for_llm(res, creq.get("ask"))
    exp = llm(EXPLAIN_SYS, "[질문]\n" + text + "\n\n[결과 요약]\n" + summ, temperature=0.2, on_token=lambda t: emit({"token": t, "ch": "explain"}))
    exp = re.sub(r"<think>.*?</think>", "", exp, flags=re.S).strip()
    stray = sorted(numbers_of(exp) - numbers_of(summ + " " + text) - {"10"})
    return {"parsed": parsed, "result": res, "id": hid, "explain": exp, "summary": summ,
            "stray_numbers": stray}


def meta():
    return {"version": VERSION, "nuclides": S.nuclides(), "materials": {k: {"name": v["name"], "rho": v["rho"], "buildup": v["buildup"]} for k, v in S.XCOM.items() if k != "air"},
            "act_units": list(S.ACT_UNITS), "len_units": list(S.LEN_UNITS), "rate_units": list(S.RATE_UNITS), "time_units": list(S.TIME_UNITS),
            "rules": S.RULES, "quant": S.QUANT, "limits": S.LIMITS, "beta": list(S.NUC["beta"]), "neutron": {k: v["note"] for k, v in S.NEUTRON.items()},
            "sources": SOURCES, "disclaimer": DISCLAIMER, "model": MODEL, "llm": llm_ok()}


# ── HTTP ───────────────────────────────────────────────────────────────
HTML = read(os.path.join(ROOT, "ui.html")) if os.path.exists(os.path.join(ROOT, "ui.html")) else "ui.html 없음"

# ── 저작권 표기 (LICENSE·NOTICE 참고) ─────────────────────────────────────
_SIG = __import__("base64").b64decode("wqkgMjAyNiBnZ2dnODY1NyDCtyBkb25nanVraW0uZGV2QGdtYWlsLmNvbQ==").decode()
_SIG_A = __import__("base64").b64decode("Z2dnZzg2NTcgPGRvbmdqdWtpbS5kZXZAZ21haWwuY29tPg==").decode()


def signed(html):
    """화면에 저작권 표기를 붙인다. ui.html 에서 지워져도 서버가 내보낼 때 다시 붙는다."""
    name, mail = _SIG.split(" · ")
    if 'name="author"' not in html:
        meta_tag = f'<meta name="author" content="{name[7:]} <{mail}>">'
        html = html.replace("<head>", "<head>" + meta_tag, 1) if "<head>" in html else meta_tag + html
    if "data-sig" not in html:
        tag = (f'<!-- {_SIG} --><div data-sig title="{mail}" style="text-align:center;font-size:11px;color:#9aa0a6;'
               f'opacity:.55;margin:28px 0 8px">{name}</div>')
        html = html.replace("</body>", tag + "</body>", 1) if "</body>" in html else html + tag
    return html


class H(BaseHTTPRequestHandler):
    def log_message(self, fmt, *a):
        if "/api/ask" in (a[0] if a else ""):
            super().log_message(fmt, *a)

    def _send(self, body, ctype="application/json", code=200, extra=None):
        b = body if isinstance(body, bytes) else json.dumps(body, ensure_ascii=False).encode()
        self.send_response(code)
        self.send_header("X-Author", _SIG_A)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(b)

    def do_GET(self):
        path = self.path.split("?")[0]
        try:
            if path == "/api/health":
                return self._send({"ok": True})
            if path == "/api/meta":
                return self._send(meta())
            if path == "/api/history":
                return self._send(hist_list())
            if path == "/api/validation":
                return self._send(validation.run())
            m = re.fullmatch(r"/api/history/([\w-]+)", path)
            if m:
                return self._send(json.loads(read(hist_path(m.group(1)))))
            m = re.fullmatch(r"/api/report/([\w-]+)\.md", path)
            if m:
                rec = json.loads(read(hist_path(m.group(1))))
                if rec["kind"] not in ("calc", "ask"):
                    raise ValueError("보고서는 감마 계산만")
                return self._send(report_md(rec["result"]).encode(), "text/markdown; charset=utf-8",
                                  extra={"Content-Disposition": f'attachment; filename="shield-{m.group(1)}.md"'})
            if path == "/api/nuclide":
                n = (self.path.split("n=", 1) + [""])[1].split("&")[0]
                g = S.NUC["gamma"].get(n)
                if not g:
                    raise FileNotFoundError
                return self._send(dict(g, gamma_const=S.gamma_constant(S.photon_lines(n, cutoff_kev=20)),
                                       gamma_const10=S.gamma_constant(S.photon_lines(n, cutoff_kev=10))))
            self._send(signed(HTML).encode(), "text/html; charset=utf-8")
        except (FileNotFoundError, ValueError):
            self._send({"error": "없음"}, code=404)
        except Exception as e:
            self._send({"error": f"{type(e).__name__}: {e}"}, code=500)

    def do_POST(self):
        try:
            req = json.loads(self.rfile.read(int(self.headers.get("Content-Length") or 0)) or b"{}")
        except ValueError:
            return self._send({"error": "잘못된 요청"}, code=400)
        try:
            if self.path == "/api/calc":
                res = calc(req)
                res["id"] = hist_save("calc", req, res, title_of(res))
                return self._send(res)
            if self.path == "/api/other":
                res = calc_other(req)
                res["id"] = hist_save(res["mode"], req, res, (req.get("nuclide") or req.get("source") or "") + (" 베타" if res["mode"] == "beta" else " 중성자"))
                return self._send(res)
            if self.path == "/api/history/delete":
                p = hist_path(req.get("id"))
                if os.path.exists(p):
                    os.remove(p)
                return self._send({"ok": True})
        except ValueError as e:
            return self._send({"error": str(e)}, code=400)
        except Exception as e:
            return self._send({"error": f"{type(e).__name__}: {e}"}, code=500)
        if self.path != "/api/ask":
            return self._send({"error": "없는 경로"}, code=404)
        self.send_response(200)
        self.send_header("X-Author", _SIG_A)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()

        def emit(ev):
            self.wfile.write(f"data: {json.dumps(ev, ensure_ascii=False)}\n\n".encode())
            self.wfile.flush()
        try:
            emit({"done": ask(req, emit)})
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as e:
            try:
                emit({"error": str(e) if isinstance(e, (ValueError, RuntimeError)) else f"{type(e).__name__}: {e}"})
            except OSError:
                pass


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--cli":
        r = calc(json.loads(sys.argv[2]))
        print(report_md(r))
        sys.exit(0)
    print(f"shield local → http://localhost:{PORT}  (llm={LLM_API} {LLM_BASE} {MODEL}, workspace={WS})  {_SIG}")
    ThreadingHTTPServer(("0.0.0.0", PORT), H).serve_forever()

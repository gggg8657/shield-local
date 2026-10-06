#!/usr/bin/env python3
"""LLM 없이 검증 (가짜 LLM): 단위 변환 → 감쇠계수 보간 → G-P 축적인자 → 점선원 선량률 → 역산 두께·HVL/TVL → 누적선량
→ 문헌 기준값 대조(validation.py, 허용오차) → 베타·중성자 → LLM 입력 검증(숫자 날조 감지) → 이력 → HTTP(API·SSE·보고서) → ui.html.
WORKSPACE 는 임시 폴더로 바꿔 실데이터 폴더에 흔적을 남기지 않는다.   python3 selftest.py"""
import json
import math
import os
import shutil
import sys
import tempfile
import threading
import urllib.request

TMP = tempfile.mkdtemp(prefix="shield-selftest-")
os.environ["WORKSPACE"] = TMP
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import app  # noqa: E402
import shield as S  # noqa: E402
import validation  # noqa: E402


def near(a, b, rel):
    return abs(a / b - 1) <= rel


def fake(system, user, temperature=0.0, on_token=lambda t: None, json_mode=False):
    if "입력 변환기" in system:
        if "에서 5 μSv" in user:  # 질문에 없는 숫자(20 mCi)를 지어냄 → 검출돼야 함
            out = '{"nuclide":"Cs137","activity":{"value":20,"unit":"mCi"},"distance":{"value":30,"unit":"cm"},"layers":[],"target":{"value":5,"unit":"μSv/h"},"solve_material":"lead","ask":"thickness","missing":[]}'
        else:
            out = ('```json\n{"nuclide":"137Cs","activity":{"value":10,"unit":"mCi"},"distance":{"value":30,"unit":"cm"},"layers":[],'
                   '"target":{"value":2.5,"unit":"μSv/h"},"solve_material":"lead","quantity":"hstar","time":null,"ask":"thickness","missing":[]}\n```')
    else:
        out = "필요한 납 두께는 4.844 cm 입니다. 비차폐 선량률은 380.2 μSv/h 이고 99.9 라는 엉뚱한 숫자도 있습니다."
    for i in range(0, len(out), 9):
        on_token(out[i:i + 9])
    return out


app.llm = fake
try:
    assert app.WS == TMP, "WORKSPACE 가 임시 폴더가 아님"

    # 1) 단위 변환
    assert S.to_bq(1, "Ci") == 3.7e10 and near(S.to_bq(10, "mCi"), 3.7e8, 1e-12) and near(S.to_bq(1, "μCi"), 3.7e4, 1e-12)
    assert S.to_bq(2, "GBq") == 2e9 and S.to_cm(5, "mm") == 0.5 and S.to_cm(1.5, "m") == 150
    assert near(S.rate_to_base(2.5, "µSv/h"), 2.5e-6, 1e-12) and S.hours(30, "min") == 0.5
    for bad in (lambda: S.to_bq(1, "Ci/kg"), lambda: S.to_cm(1, "ft")):
        try:
            bad()
            raise AssertionError("잘못된 단위를 받아들임")
        except ValueError:
            pass

    # 2) 감쇠계수 보간 (NIST 격자점·보간점·흡수끝)
    assert near(S.mu_rho("lead", 1.0, True), 0.07102, 1e-4) and near(S.mu_rho("water", 1.0, True), 0.07072, 1e-3)
    assert near(S.mu_rho("lead", 0.6617, True), 0.1102, 3e-3) and near(S.mu_rho("lead", 0.6617), 0.1035, 3e-3)
    assert S.mu_rho("lead", 0.6617) < S.mu_rho("lead", 0.6617, True)  # 간섭성 제외 값이 더 작다
    below, above = S.mu_rho("lead", 0.0879, True), S.mu_rho("lead", 0.0881, True)
    assert below < 2.0 and above > 7.0, (below, above)  # 납 K 흡수끝 88.0 keV 를 건너 보간하지 않음
    assert near(S.mu_rho("lead", 0.088, True), 1.91, 0.01)  # 끝 에너지 정확히 → 끝 아래 값(보수적)
    assert near(S.muen_air(1.25), 0.02666, 1e-3) and S.h_per_ka(0.005) == 0 and near(S.h_per_ka(0.6617), 1.203, 0.01)

    # 3) G-P 축적인자: Table 3 재구성, 경계
    for m, e, x, b in [("lead", 1.0, 4, 2.10), ("concrete", 1.0, 4, 6.42), ("water", 2.0, 7, 8.65), ("concrete", 8.0, 20, 8.31)]:
        assert near(S.buildup(m, e, x), b, 0.03), (m, e, x, S.buildup(m, e, x), b)
    assert S.buildup("lead", 1.0, 0) == 1.0 and S.buildup("lead", 1.0, 100) == S.buildup("lead", 1.0, 40)
    assert S.buildup("concrete", 1.0, 2) > S.buildup("lead", 1.0, 2)  # 저Z 가 산란 축적이 크다

    # 4) 점선원 선량률: Co-60 1 Ci 1 m, 역제곱, 차폐 단조
    co = S.photon_lines("Co-60")
    k = S.point_rate(co, 3.7e10, 100, [])["kerma"] * 1e3
    assert near(k, 11.43, 0.03), k  # Γ 309 μGy·m²/(GBq·h) × 37 GBq (Ninković 2012)
    assert near(S.point_rate(co, 1e9, 200, [])["kerma"] * 4, S.point_rate(co, 1e9, 100, [])["kerma"], 1e-12)
    pb = lambda t: S.point_rate(co, 1e9, 100, S.norm_layers([{"material": "lead", "thickness": t, "unit": "cm"}]))["hstar"]
    assert pb(0) > pb(1) > pb(5) > pb(10) > 0
    two = S.norm_layers([{"material": "concrete", "thickness": 20, "unit": "cm"}, {"material": "lead", "thickness": 2, "unit": "cm"}])
    rmax, rlast, rnone = (S.point_rate(co, 1e9, 100, two, r)["kerma"] for r in ("max", "last", "none"))
    assert rmax >= rlast >= rnone, (rmax, rlast, rnone)
    assert len(S.photon_lines("Am-241", cutoff_kev=10)) > len(S.photon_lines("Am-241", cutoff_kev=20)) > 0

    # 5) 역산: 목표 선량률 → 두께 (되대입), HVL/TVL
    cs = S.photon_lines("Cs-137")
    a = S.to_bq(10, "mCi")
    t = S.solve_thickness(cs, a, 30, [], "lead", None, 2.5e-6, "hstar")
    back = S.point_rate(cs, a, 30, [{"mat": "lead", "t_cm": t, "rho": 11.35}])["hstar"]
    assert near(back, 2.5e-6, 1e-4), back
    assert 3.0 < t < 6.0, t
    assert S.solve_thickness(cs, a, 30, [], "lead", None, 1.0, "hstar") == 0.0  # 이미 목표 이하
    hv = S.hvl_tvl(cs, "lead", quantity="kerma")
    assert hv["hvl_narrow"] < hv["hvl_broad"] and near(hv["hvl_narrow"], math.log(2) / (S.mu_rho("lead", 0.6617) * 11.35), 0.05)

    # 6) 누적선량 (붕괴)
    assert S.cumulative(1.0, 2.0, None) == 2.0
    f18 = S.half_life_s("F-18")
    assert near(S.cumulative(1.0, 1e6, f18), f18 / 3600 / math.log(2), 1e-6)

    # 7) 문헌 기준값 대조 — 모두 허용오차 이내
    rows = validation.run()
    bad = [r for r in rows if not r["ok"]]
    assert not bad, bad
    assert len(rows) >= 40

    # 8) 베타·중성자
    b = S.beta_info("Sr-90")
    y90 = [x for x in b["branches"] if x["origin"] == "Y-90"][0]
    assert near(y90["Emax_MeV"], 2.2787, 1e-3) and 1.0 < y90["range_gcm2"] < 1.2  # Katz–Penfold ≈ 1.1 g/cm²
    n = S.neutron_rate("Cf-252", 1, "ug", 100)
    assert near(n["S_n_per_s"], 2.31e6, 0.01) and near(n["hstar"] * 1e6, 25.6, 0.02), n

    # 9) LLM 입력 검증: 핵종 정규화·μ 단위·지어낸 숫자 감지
    p = app.validate_parsed({"nuclide": "137Cs", "activity": {"value": 10, "unit": "mCi"}, "distance": {"value": 30, "unit": "cm"},
                             "target": {"value": 2.5, "unit": "μSv/h"}}, "Cs-137 10 mCi 30 cm 2.5 μSv/h 납 몇 mm?")
    assert p["request"]["source"]["nuclide"] == "Cs-137" and p["request"]["target"]["unit"] == "uSv/h" and not p["issues"], p
    p = app.validate_parsed({"nuclide": "Cs-137", "activity": {"value": 20, "unit": "mCi"}, "distance": {"value": 30, "unit": "cm"}}, "Cs-137 10 mCi 30 cm")
    assert any("질문에 없음" in i for i in p["issues"]), p
    p = app.validate_parsed({"nuclide": "Xx-999", "activity": None, "distance": None}, "아무거나")
    assert any("지원하지 않는" in i for i in p["issues"]) and any("거리" in i for i in p["issues"])

    # 10) calc 전체 + 이력 + 보고서
    req = {"source": {"nuclide": "Cs-137", "activity": {"value": 10, "unit": "mCi"}}, "distance": {"value": 30, "unit": "cm"},
           "layers": [{"material": "lead", "thickness": 10, "unit": "mm"}], "target": {"value": 2.5, "unit": "uSv/h"},
           "solve_material": "lead", "time": {"value": 2, "unit": "h"}}
    res = app.calc(req)
    assert res["steps"] and res["curve_t"]["points"] and res["curve_r"]["points"] and res["hvl"]["tvl_eq"]
    assert res["time"]["limits"][2]["frac"] > 0 and res["disclaimer"].startswith("간이 계산")
    try:
        app.calc(dict(req, distance={"value": 5, "unit": "mm"}))
        raise AssertionError("차폐 두께 > 거리 를 받아들임")
    except ValueError:
        pass
    md = app.report_md(res)
    assert "간이 계산" in md and "데이터 출처" in md and "NIST XCOM" in md and "ANS-6.4.3" in md

    # 11) HTTP
    srv = app.ThreadingHTTPServer(("127.0.0.1", 0), app.H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"

    def post(path, body):
        return urllib.request.urlopen(urllib.request.Request(base + path, json.dumps(body).encode(), {"Content-Type": "application/json"}), timeout=60)
    html = urllib.request.urlopen(base + "/").read().decode()
    assert "data-sig" in html and 'name="author"' in html and "cdn" not in html.lower()
    r = urllib.request.urlopen(base + "/api/health")
    assert r.headers["X-Author"]
    mt = json.load(urllib.request.urlopen(base + "/api/meta"))
    assert "Cs-137" in mt["nuclides"] and "lead" in mt["materials"] and mt["disclaimer"]
    j = json.load(post("/api/calc", req))
    assert j["id"] and near(j["shielded"]["hstar"], res["shielded"]["hstar"], 1e-12)
    md = urllib.request.urlopen(base + f"/api/report/{j['id']}.md").read().decode()
    assert md.startswith("# 방사선 차폐")
    h = json.load(urllib.request.urlopen(base + "/api/history"))
    assert h and h[0]["id"] == j["id"]
    assert json.load(urllib.request.urlopen(base + "/api/history/" + j["id"]))["request"]["source"]["nuclide"] == "Cs-137"
    assert all(r["ok"] for r in json.load(urllib.request.urlopen(base + "/api/validation")))
    assert json.load(urllib.request.urlopen(base + "/api/nuclide?n=Co-60"))["gamma_const"] > 300
    assert json.load(post("/api/other", {"mode": "beta", "nuclide": "P-32"}))["branches"][0]["Emax_MeV"] > 1.7
    try:
        post("/api/calc", {"source": {"nuclide": "Cs-137"}})
        raise AssertionError("방사능 없이 계산됨")
    except urllib.error.HTTPError as e:
        assert e.code == 400

    def sse(body):
        evs = []
        for line in post("/api/ask", body):
            line = line.decode().strip()
            if line.startswith("data:"):
                evs.append(json.loads(line[5:]))
        return evs
    evs = sse({"question": "Cs-137 10 mCi 를 30 cm 에서 2.5 μSv/h 이하로 하려면 납 몇 mm?"})
    done = [e for e in evs if "done" in e][0]["done"]
    assert any("token" in e for e in evs) and done["result"]["target"]["t_cm"] > 3
    assert "99.9" in done["stray_numbers"] and "4.844" not in done["stray_numbers"], done["stray_numbers"]
    evs = sse({"question": "Cs-137 10 mCi 30 cm 에서 5 μSv/h 되려면?"})
    done = [e for e in evs if "done" in e][0]["done"]
    assert any("질문에 없음" in i for i in done["parsed"]["issues"])
    post("/api/history/delete", {"id": j["id"]})
    assert j["id"] not in [x["id"] for x in json.load(urllib.request.urlopen(base + "/api/history"))]
    srv.shutdown()
    print(f"selftest OK — 단위·보간·축적인자·선량률·역산·누적·문헌대조 {len(rows)}건·베타·중성자·LLM 검증·HTTP")
finally:
    shutil.rmtree(TMP, ignore_errors=True)

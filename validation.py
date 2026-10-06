"""문헌 기준값 대조 — selftest·화면('데이터·검증' 탭)·README 가 같은 표를 쓴다.

  python3 validation.py        # 대조표를 Markdown 으로 출력

tol 은 '이 범위를 넘으면 selftest 실패' 기준. 문헌끼리의 차이·데이터 판(版) 차이를 감안해 정했다.
"""
import shield as S

NIST = "NIST XCOM / Hubbell-Seltzer 표"
NINK = "Ninković & Adrović (2012), δ=20 keV"
RHH = "Radiological Health Handbook Γ(R·cm²/mCi·h) × 23.69"  # 1 R = 8.764 mGy 공기커마, 1 mCi = 0.037 GBq
NCRP = "NCRP 49 (1976) 넓은 빔 평형 TVL"
TR = "NUREG/CR-5740 Table 3"


def _mu(mat, e, coh=True):
    return lambda: S.mu_rho(mat, e, coh)


def _gam(nuc, cut=20.0):
    return lambda: S.gamma_constant(S.photon_lines(nuc, cutoff_kev=cut))


def _B(mat, e, x):
    return lambda: S.buildup(mat, e, x)


def _tvl_eq(nuc, mat):
    def f():
        lines = S.photon_lines(nuc)
        b = S.point_rate(lines, 1.0, 100.0, [], "max")["kerma"]
        t2, t3 = (S.solve_thickness(lines, 1.0, 100.0, [], mat, None, b * k, "kerma", "max") for k in (1e-2, 1e-3))
        return t3 - t2
    return f


def _co60_1ci():
    return S.point_rate(S.photon_lines("Co-60"), S.to_bq(1, "Ci"), 100.0, [])["kerma"] * 1e3


# (구분, 항목, 계산 함수, 문헌값, 단위, 출처, 허용 상대오차, 메모)
CASES = [
    ("감쇠계수", "납 1 MeV μ/ρ", _mu("lead", 1.0), 0.07102, "cm²/g", NIST, 0.005, "간섭성 포함(표 값)"),
    ("감쇠계수", "납 0.5 MeV μ/ρ", _mu("lead", 0.5), 0.1614, "cm²/g", NIST, 0.005, ""),
    ("감쇠계수", "납 661.7 keV μ/ρ (보간)", _mu("lead", 0.6617), 0.1102, "cm²/g", "NIST XCOM 직접 조회", 0.005, "격자 보간 검사"),
    ("감쇠계수", "납 661.7 keV μ/ρ 간섭성 제외 (보간)", _mu("lead", 0.6617, False), 0.1035, "cm²/g", "NIST XCOM 직접 조회", 0.005, "차폐 계산에 쓰는 값"),
    ("감쇠계수", "물 1 MeV μ/ρ", _mu("water", 1.0), 0.07072, "cm²/g", NIST, 0.005, ""),
    ("감쇠계수", "보통 콘크리트 1 MeV μ/ρ", _mu("concrete", 1.0), 0.06495, "cm²/g", NIST, 0.005, "NIST Concrete, Ordinary"),
    ("감쇠계수", "철 1 MeV μ/ρ", _mu("iron", 1.0), 0.05995, "cm²/g", NIST, 0.005, ""),
    ("감쇠계수", "텅스텐 1 MeV μ/ρ", _mu("tungsten", 1.0), 0.06618, "cm²/g", NIST, 0.005, ""),
    ("감쇠계수", "공기 μen/ρ 1.25 MeV", lambda: S.muen_air(1.25), 0.02666, "cm²/g", NIST, 0.005, ""),
    ("축적인자", "납 1 MeV, 4 mfp", _B("lead", 1.0, 4), 2.10, "", TR, 0.03, "G-P 재구성"),
    ("축적인자", "납 1 MeV, 10 mfp", _B("lead", 1.0, 10), 3.37, "", TR, 0.03, ""),
    ("축적인자", "납 0.5 MeV, 4 mfp", _B("lead", 0.5, 4), 1.53, "", TR, 0.03, ""),
    ("축적인자", "콘크리트 1 MeV, 4 mfp", _B("concrete", 1.0, 4), 6.42, "", TR, 0.03, ""),
    ("축적인자", "콘크리트 0.5 MeV, 10 mfp", _B("concrete", 0.5, 10), 36.4, "", TR, 0.03, ""),
    ("축적인자", "콘크리트 8 MeV, 20 mfp", _B("concrete", 8.0, 20), 8.31, "", TR, 0.03, "계수 전사 오류 정정 확인용"),
    ("축적인자", "물 1.5 MeV, 10 mfp", _B("water", 1.5, 10), 16.7, "", TR, 0.03, ""),
    ("축적인자", "철 2 MeV, 10 mfp", _B("iron", 2.0, 10), 10.8, "", TR, 0.03, ""),
    ("선량률", "Co-60 1 Ci, 1 m 비차폐 공기커마율", _co60_1ci, 309.0 * 37 / 1000, "mGy/h", NINK + " Γ×37 GBq", 0.03, "H*(10) 로는 약 13.1 mSv/h"),
    ("선량률", "Co-60 1 Ci, 1 m 비차폐 공기커마율 ", _co60_1ci, 1.32 * 8.764, "mGy/h", "RHH 1.32 R/h × 8.764 mGy/R", 0.03, ""),
    ("Γ 상수", "Co-60", _gam("Co-60"), 309.0, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Co-60 ", _gam("Co-60"), 13.2 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "Cs-137 (Ba-137m)", _gam("Cs-137"), 3.3 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "Cs-137 (Ba-137m) ", _gam("Cs-137"), 82.10, "μGy·m²/(GBq·h)", NINK, 0.08,
     "−6 %: 문헌이 661.7 keV 방출률을 더 크게 쓴 것으로 보임(DDEP 85.01 %). RHH 와는 −2 %"),
    ("Γ 상수", "Ir-192", _gam("Ir-192"), 109.1, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Ir-192 ", _gam("Ir-192"), 4.8 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "Na-22", _gam("Na-22"), 12.0 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "Mn-54", _gam("Mn-54"), 4.7 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "I-131", _gam("I-131"), 52.20, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "I-131 ", _gam("I-131"), 2.2 * 23.69, "μGy·m²/(GBq·h)", RHH, 0.05, ""),
    ("Γ 상수", "F-18", _gam("F-18"), 135.1, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Se-75", _gam("Se-75"), 48.25, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Eu-152", _gam("Eu-152"), 148.9, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Eu-154", _gam("Eu-154"), 159.2, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Co-58", _gam("Co-58"), 129.0, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Ga-68", _gam("Ga-68"), 129.0, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Fe-59", _gam("Fe-59"), 145.9, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Cr-51", _gam("Cr-51"), 4.22, "μGy·m²/(GBq·h)", NINK, 0.03, ""),
    ("Γ 상수", "Tc-99m", _gam("Tc-99m"), 14.10, "μGy·m²/(GBq·h)", NINK, 0.06, "저에너지 X선 비중 큼"),
    ("Γ 상수", "Tl-201", _gam("Tl-201"), 10.22, "μGy·m²/(GBq·h)", NINK, 0.06, "저에너지 X선 비중 큼"),
    ("Γ 상수", "Co-57", _gam("Co-57"), 14.11, "μGy·m²/(GBq·h)", NINK, 0.08, "저에너지 X선 비중 큼"),
    ("Γ 상수", "I-125", _gam("I-125"), 37.73, "μGy·m²/(GBq·h)", NINK, 0.10, "27~35 keV — X선 데이터 판 차이"),
    ("Γ 상수", "Am-241 (δ=22 keV)", _gam("Am-241", 22.0), 3.97, "μGy·m²/(GBq·h)", NINK, 0.08,
     "δ=20 keV 면 Np Lγ X선(21.2 keV, 4.8 %)이 들어가 5.8 — 문헌은 이를 뺀 값으로 보임"),
    ("평형 TVL", "Cs-137 · 납", _tvl_eq("Cs-137", "lead"), 2.1, "cm", NCRP, 0.10, "1/100→1/1000 구간 두께"),
    ("평형 TVL", "Cs-137 · 콘크리트", _tvl_eq("Cs-137", "concrete"), 15.7, "cm", NCRP, 0.10, ""),
    ("평형 TVL", "Co-60 · 납", _tvl_eq("Co-60", "lead"), 4.0, "cm", NCRP, 0.10, ""),
    ("평형 TVL", "Co-60 · 콘크리트", _tvl_eq("Co-60", "concrete"), 20.6, "cm", NCRP, 0.10, ""),
    ("평형 TVL", "Ir-192 · 납", _tvl_eq("Ir-192", "lead"), 1.9, "cm", NCRP, 0.12, ""),
    ("평형 TVL", "Ir-192 · 콘크리트", _tvl_eq("Ir-192", "concrete"), 14.7, "cm", NCRP, 0.12, "NCRP 값은 밀도 2.35 g/cm³"),
]


def run():
    out = []
    for grp, name, fn, ref, unit, src, tol, note in CASES:
        v = fn()
        err = v / ref - 1
        out.append({"group": grp, "name": name.strip(), "calc": v, "ref": ref, "unit": unit, "source": src, "tol": tol,
                    "err": err, "ok": abs(err) <= tol, "note": note})
    return out


def markdown(rows=None):
    rows = rows or run()
    L = ["| 구분 | 항목 | 계산값 | 문헌값 | 단위 | 오차 | 허용 | 출처 | 메모 |", "|---|---|---:|---:|---|---:|---:|---|---|"]
    for r in rows:
        L.append(f"| {r['group']} | {r['name']} | {r['calc']:.4g} | {r['ref']:.4g} | {r['unit']} | {100 * r['err']:+.1f}% | ±{100 * r['tol']:g}% | {r['source']} | {r['note']} |")
    return "\n".join(L)


if __name__ == "__main__":
    rows = run()
    print(markdown(rows))
    print(f"\n{sum(r['ok'] for r in rows)}/{len(rows)} 허용오차 이내")

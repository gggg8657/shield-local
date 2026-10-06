#!/usr/bin/env python3
"""data/*.json 을 공개 원본에서 다시 만든다 (인터넷 필요 — 폐쇄망 서버에서는 실행할 필요 없음, 결과 JSON 을 동봉).

  python3 tools/build_data.py            # data/xcom.json, data/air_muen.json, data/nuclides.json 재생성

원본
  - 감쇠계수: NIST XCOM (Berger, Hubbell, Seltzer et al., NIST Standard Reference Database 8, XCOM v3.1)
      https://physics.nist.gov/cgi-bin/Xcom/xcom3_1  — 원소별로 받아 혼합물은 질량분율 가중합(XCOM 과 같은 방식)
      혼합물 조성: NIST X-Ray Mass Attenuation Coefficients Table 2 (Hubbell & Seltzer, NISTIR 5632)
  - 공기 질량에너지흡수계수 μen/ρ: NIST X-Ray Mass Attenuation Coefficients (Hubbell & Seltzer, NISTIR 5632), Air, Dry
      https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
  - 붕괴 데이터: DDEP 권고값 (LNHB/CEA, Laboratoire National Henri Becquerel) LARA 파일
      http://www.lnhb.fr/nuclides/<핵종>.lara.txt
G-P 축적인자(data/buildup_ans643.json)와 ICRP 74 환산계수는 표를 직접 대조해 손으로 넣었다(이 스크립트 대상 아님).
"""
import json
import math
import os
import re
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data")
UA = {"User-Agent": "Mozilla/5.0 (shield-local build_data)"}


def get(url, data=None):
    req = urllib.request.Request(url, data, UA)
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", "replace")


# ── 재질 ─────────────────────────────────────────────────────────────────
# 원소 재질: Z.  혼합물: NIST Table 2 질량분율.  rho = 기본 밀도(g/cm³, 화면에서 바꿀 수 있음)
MATERIALS = {
    "lead": {"name": "납 (Pb)", "Z": 82, "rho": 11.35, "buildup": "lead"},
    "tungsten": {"name": "텅스텐 (W)", "Z": 74, "rho": 19.30, "buildup": "tungsten"},
    "iron": {"name": "철 (Fe)", "Z": 26, "rho": 7.874, "buildup": "iron"},
    "ss304": {"name": "스테인리스강 SS304", "rho": 8.00, "buildup": "iron",
              "frac": {24: 0.19, 25: 0.02, 26: 0.695, 28: 0.095}},
    "copper": {"name": "구리 (Cu)", "Z": 29, "rho": 8.96, "buildup": "copper"},
    "aluminium": {"name": "알루미늄 (Al)", "Z": 13, "rho": 2.699, "buildup": "aluminium"},
    "concrete": {"name": "보통 콘크리트", "rho": 2.30, "buildup": "concrete",
                 "frac": {1: 0.022100, 6: 0.002484, 8: 0.574930, 11: 0.015208, 12: 0.001266, 13: 0.019953,
                          14: 0.304627, 19: 0.010045, 20: 0.042951, 26: 0.006435}},
    "water": {"name": "물", "rho": 1.00, "buildup": "water", "frac": {1: 0.111898, 8: 0.888102}},
    "polyethylene": {"name": "폴리에틸렌 (PE)", "rho": 0.94, "buildup": "water", "frac": {1: 0.143716, 6: 0.856284}},
    "air": {"name": "공기 (건조)", "rho": 0.001205, "buildup": "air",
            "frac": {6: 0.000124, 7: 0.755268, 8: 0.231781, 18: 0.012827}},
}

# 조밀 격자: 10 keV ~ 20 MeV, 10년당 40점 (log 간격 5.9 %) — log-log 보간 오차 ≪ 0.5 %
DENSE = sorted({round(10 ** (k / 40) * 0.01, 6) for k in range(0, int(40 * math.log10(2000)) + 1)} | {0.01, 20.0})
ROW = re.compile(r"^\s*(?:(\d+)\s*(?:&#160;|\s)\s*([KLM]\s*\d?)\s+)?(\d\.\d{3}E[+-]\d\d)((?:\s+\d\.\d{3}E[+-]\d\d){7})\s*$")


def xcom_element(z):
    """[(E, mu_with_coh, mu_without_coh, edge_label)] 0.01~20 MeV, 표준 격자(흡수끝 포함) + 조밀 격자"""
    rows = []
    for k in range(0, len(DENSE), 40):  # XCOM 은 한 번에 받는 추가 에너지 수가 제한됨 — 나눠 받고 표준 격자(흡수끝 포함)는 첫 번에만
        q = {"ZNum": str(z), "OutOpt": "PIC", "NumAdd": "1", "Energies": "\r\n".join(f"{e:g}" for e in DENSE[k:k + 40]),
             "WindowXmin": "0.001", "WindowXmax": "100000", "ResizeFlag": "on"} | ({"Output": "on"} if k == 0 else {})
        txt = get("https://physics.nist.gov/cgi-bin/Xcom/xcom3_1", urllib.parse.urlencode(q).encode())
        if "XCOM: Error" in txt:
            raise RuntimeError(f"XCOM 오류 Z={z}")
        for line in re.sub(r"<[^>]*>", " ", txt).split("\n"):
            m = ROW.match(line.replace("\xa0", " "))
            if not m:
                continue
            e = float(m.group(3))
            v = [float(x) for x in m.group(4).split()]
            if 0.01 <= e <= 20.0:
                rows.append((e, v[5], v[6], (m.group(2) or "").replace(" ", "")))
    rows.sort(key=lambda r: (r[0], bool(r[3])))  # 같은 에너지면 흡수끝 아래 값(라벨 없음) → 위 값(라벨) 순
    # 같은 에너지 중복(표준 격자와 조밀 격자 겹침) 제거 — 흡수끝(라벨 있는 줄)은 남김
    out = []
    for r in rows:
        if out and abs(out[-1][0] - r[0]) < 1e-9 and not r[3] and not out[-1][3]:
            continue
        out.append(r)
    return out


def interp_edge(rows, e, col):
    """원소 표에서 e 의 값 (log-log). e 가 흡수끝이면 끝 아래 값."""
    for i in range(len(rows) - 1):
        e0, e1 = rows[i][0], rows[i + 1][0]
        if e0 <= e <= e1 and e1 > e0:
            y0, y1 = rows[i][col], rows[i + 1][col]
            return math.exp(math.log(y0) + (math.log(y1) - math.log(y0)) * (math.log(e) - math.log(e0)) / (math.log(e1) - math.log(e0)))
        if abs(e - e0) < 1e-12:
            return rows[i][col]
    return rows[-1][col]


def build_xcom():
    zs = sorted({m["Z"] for m in MATERIALS.values() if "Z" in m} | {z for m in MATERIALS.values() for z in m.get("frac", {})})
    elem = {}
    for z in zs:
        elem[z] = xcom_element(z)
        print(f"  XCOM Z={z}: {len(elem[z])} 행")
    mats = {}
    for key, m in MATERIALS.items():
        if "Z" in m:
            rows = [[r[0], r[1], r[2]] + ([r[3]] if r[3] else []) for r in elem[m["Z"]]]
        else:
            grid = sorted({r[0] for z in m["frac"] for r in elem[z] if not r[3]})
            rows = [[e, round(sum(w * interp_edge(elem[z], e, 1) for z, w in m["frac"].items()), 6),
                     round(sum(w * interp_edge(elem[z], e, 2) for z, w in m["frac"].items()), 6)] for e in grid]
        mats[key] = {k: v for k, v in m.items() if k != "frac"} | ({"composition": {str(z): w for z, w in m["frac"].items()}} if "frac" in m else {}) | {"table": rows}
    doc = {"_source": "NIST XCOM v3.1 (Berger, Hubbell, Seltzer, Chang, Coursey, Sukumar, Zucker, Olsen; NIST SRD 8), "
                      "https://physics.nist.gov/cgi-bin/Xcom/xcom3_1 ; 혼합물 = 원소 μ/ρ 의 질량분율 가중합, 조성은 NIST Table 2 (NISTIR 5632). "
                      "SS304 조성은 공칭값(Fe 69.5, Cr 19, Ni 9.5, Mn 2 wt%).",
           "_columns": ["E (MeV)", "mu/rho with coherent (cm2/g)", "mu/rho without coherent (cm2/g)", "흡수끝 라벨(있으면 끝 위 값)"],
           "materials": mats}
    json.dump(doc, open(os.path.join(OUT, "xcom.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=0)


def build_air():
    txt = get("https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html")
    pre = re.search(r"<PRE>(.*?)</PRE>", txt, re.S | re.I).group(1)
    rows = []
    for line in re.sub(r"<[^>]*>", " ", pre).split("\n"):
        m = re.match(r"^\s*(?:[KLM]\d?\s+)?(\d\.\d+E[+-]\d\d)\s+(\d\.\d+E[+-]\d\d)\s+(\d\.\d+E[+-]\d\d)\s*$", line)
        if m and 0.01 <= float(m.group(1)) <= 20:
            rows.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
    json.dump({"_source": "NIST X-Ray Mass Attenuation Coefficients, Table 4 Air, Dry (near sea level) — J. H. Hubbell & S. M. Seltzer, NISTIR 5632 (1995, 2004 web ed.)",
               "_columns": ["E (MeV)", "mu/rho (cm2/g)", "mu_en/rho (cm2/g)"], "table": rows},
              open(os.path.join(OUT, "air_muen.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=0)
    print(f"  air μen/ρ: {len(rows)} 행")


# ── 핵종 ─────────────────────────────────────────────────────────────────
# (키, LARA 파일들 — 평형 딸핵종은 [파일, 분기비] 로 더함, 화면 이름, 메모)
NUCLIDES = [
    ("Co-60", [("Co-60", 1)], ""), ("Cs-137", [("Cs-137", 1)], "Ba-137m 661.7 keV 포함(Cs-137 붕괴당)"),
    ("Ir-192", [("Ir-192", 1)], ""), ("Am-241", [("Am-241", 1)], "Np L X선은 13~21 keV — 차단에너지 설정에 주의"),
    ("Ba-133", [("Ba-133", 1)], ""), ("Na-22", [("Na-22", 1)], "511 keV 소멸광자 180.7 % 포함"),
    ("Mn-54", [("Mn-54", 1)], ""), ("Eu-152", [("Eu-152", 1)], ""), ("Eu-154", [("Eu-154", 1)], ""),
    ("I-131", [("I-131", 1)], ""), ("I-125", [("I-125", 1)], "광자 에너지 27~35 keV"),
    ("Tc-99m", [("Tc-99m", 1)], ""), ("Mo-99", [("Mo-99", 1)], "DDEP Mo-99 목록 그대로 — 140.5 keV(89.6 %, Tc-99m 평형) 포함"),
    ("Se-75", [("Se-75", 1)], ""), ("F-18", [("F-18", 1)], "511 keV 소멸광자"), ("Ga-68", [("Ga-68", 1)], ""),
    ("Ge-68", [("Ge-68", 1), ("Ga-68", 1)], "Ga-68 영속평형 포함"), ("Co-57", [("Co-57", 1)], ""), ("Co-58", [("Co-58", 1)], ""),
    ("Cd-109", [("Cd-109", 1)], "88 keV + Ag K X선"), ("Zn-65", [("Zn-65", 1)], ""), ("Cr-51", [("Cr-51", 1)], ""),
    ("Fe-59", [("Fe-59", 1)], ""), ("Lu-177", [("Lu-177", 1)], ""), ("In-111", [("In-111", 1)], ""), ("Tl-201", [("Tl-201", 1)], ""),
    ("Sm-153", [("Sm-153", 1)], ""), ("Sb-124", [("Sb-124", 1)], ""),
    ("Ra-226", [("Ra-226", 1), ("Rn-222", 1), ("Pb-214", 1), ("Bi-214", 1)], "딸핵종(Rn-222·Pb-214·Bi-214) 평형 포함 — 밀봉선원 가정"),
]
BETA = {  # 순수(또는 주) 베타 방출체: Q- (keV, LARA) ≈ 주 분기 최대에너지
    "Sr-90": ["Sr-90", "Y-90"], "Y-90": ["Y-90"], "P-32": ["P-32"], "H-3": ["H-3"], "C-14": ["C-14"], "S-35": ["S-35"],
    "Cl-36": ["Cl-36"], "Ni-63": ["Ni-63"], "Tl-204": ["Tl-204"], "Kr-85": ["Kr-85"],
}


def lara(name):
    txt = get(f"http://www.lnhb.fr/nuclides/{name}.lara.txt")
    head = {}
    lines = []
    for line in txt.splitlines():
        p = [s.strip() for s in line.split(";")]
        if len(p) >= 5 and re.fullmatch(r"[\d.]+", p[0] or "x") and p[2]:
            lines.append({"E_keV": float(p[0]), "I": float(p[2]) / 100.0, "type": p[4], "origin": p[5] if len(p) > 5 else ""})
        elif len(p) >= 2 and p[0]:
            head[p[0]] = p[1:]
    return head, lines


def build_nuclides():
    nucs = {}
    for key, parts, note in NUCLIDES:
        photons, refs, half = [], [], None
        for i, (f, br) in enumerate(parts):
            head, lines = lara(f)
            refs.append(f"{f}: {head.get('Reference', ['?'])[0]}")
            if i == 0:
                half = float(head["Half-life (s)"][0])
            for l in lines:
                if (l["type"].startswith("g") or l["type"].startswith("X")) and l["E_keV"] >= 10.0 and l["I"] * br >= 1e-5:
                    photons.append([round(l["E_keV"] / 1000, 7), round(l["I"] * br, 7), "X" if l["type"].startswith("X") else "g", f])
        photons.sort()
        nucs[key] = {"half_life_s": half, "photons": photons, "note": note, "ddep_ref": refs}
        print(f"  {key}: {len(photons)} 광자선, T½={half:.4g} s")
    betas = {}
    for key, files in BETA.items():
        qs, half = [], None
        for i, f in enumerate(files):
            head, _ = lara(f)
            q = head.get("Q-")
            if q:
                qs.append([f, float(q[0])])
            if i == 0:
                half = float(head["Half-life (s)"][0])
        betas[key] = {"half_life_s": half, "Q_keV": qs}
    json.dump({"_source": "DDEP 권고 붕괴 데이터 (Decay Data Evaluation Project; LNHB/CEA LARA 파일 http://www.lnhb.fr/nuclides/). "
                          "광자: γ·X선, 10 keV 이상, 붕괴당 방출률 1e-5 이상. 평형 딸핵종은 분기비를 곱해 합산.",
               "_columns": ["E (MeV)", "붕괴당 방출률", "g=γ / X=X선", "출처 핵종"], "gamma": nucs, "beta": betas},
              open(os.path.join(OUT, "nuclides.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=0)


if __name__ == "__main__":
    import sys
    os.makedirs(OUT, exist_ok=True)
    what = sys.argv[1:] or ["air", "xcom", "nuclides"]  # 일부만: python3 tools/build_data.py nuclides
    for w in what:
        {"air": build_air, "xcom": build_xcom, "nuclides": build_nuclides}[w]()
    print("완료 →", OUT)

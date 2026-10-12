# THIRD_PARTY_NOTICES — shield-local

이 저장소의 코드(`app.py`·`shield.py`·`validation.py`·`ui.html`·`selftest.py`·`setup.sh`·`tools/build_data.py`)는 DongJu Kim(gggg8657)이 작성했으며 MIT 입니다. Python 표준 라이브러리만 쓰며 제3자 코드·글꼴·아이콘은 동봉하지 않습니다. 아래는 `data/` 와 코드 안에 들어 있는 **제3자 수치 데이터**와 외부 서비스입니다.

조사 기준: 2026-10-12. 실제 파일(`data/*.json` 의 `_source` 필드, `tools/build_data.py`, `shield.py`·`validation.py` 상수)을 보고 각 원출처의 이용 조건을 확인했습니다.

## 구성요소

| 이름 | 버전 | 쓰임 | 라이선스·이용 조건 | 링크 | 상업 이용 | 비고 |
|---|---|---|---|---|---|---|
| NIST XCOM (Standard Reference Database 8) — 원소별 광자 질량감쇠계수 | XCOM v3.1 | **동봉 데이터** `data/xcom.json` (`tools/build_data.py` 가 받아 가공) | NIST SRD — Standard Reference Data Act 에 따라 미 상무부 장관 명의 저작권, "사전 허락 없이 복제 불가" 문구 | https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html , https://www.nist.gov/srd/public-law | **확인 필요** | [#3](https://github.com/gggg8657/shield-local/issues/3). 혼합물은 원소값의 질량분율 가중합(조성: NISTIR 5632 Table 2) |
| NIST X-Ray Mass Attenuation Coefficients (SRD 126, Hubbell & Seltzer NISTIR 5632) — 공기 μen/ρ | 1995 / 2004 웹판 | **동봉 데이터** `data/air_muen.json` | 위와 같음(NIST SRD) | https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html | **확인 필요** | [#3](https://github.com/gggg8657/shield-local/issues/3) |
| G-P 노출 축적인자 계수 (ANSI/ANS-6.4.3-1991 데이터) — D. K. Trubey, NUREG/CR-5740 (ORNL/RSIC-49/R1, 1991) | — | **동봉 데이터** `data/buildup_ans643.json` | 미국 정부 보고서(public domain). ANS 표준 문서 자체(ANS 저작물)는 쓰지 않음 | https://doi.org/10.2172/7032799 | 가능 / 전사 경로 **확인 필요** | 기계 판독 형태는 [BoykoNeov/Radiographer](https://github.com/BoykoNeov/Radiographer) 전사본(코드는 BNCL-1.0 비상업, 데이터는 원출처 조건이라고 그 NOTICE 에 적힘)에서 출발해 5 재질을 스캔과 대조. [#2](https://github.com/gggg8657/shield-local/issues/2) |
| DDEP 권고 붕괴 데이터 (LNHB/CEA LARA 파일) | 받은 날짜는 `_source` 참고 | **동봉 데이터** `data/nuclides.json` | 이용 조건 문서 **미확인** | http://www.lnhb.fr/nuclides/ | **미확인** | [#3](https://github.com/gggg8657/shield-local/issues/3). 출처 표기 "DDEP / LNHB-CEA" |
| ICRP Publication 74 (1996) Table A.21 — 공기커마 → H*(10) 환산계수 25개 | — | 코드 상수 `shield.py` `ICRP74_HK` | ICRP 저작물(Annals of the ICRP). 표 인용 조건 **미확인** | https://www.icrp.org/ | **확인 필요** | [#1](https://github.com/gggg8657/shield-local/issues/1). "공기커마율" 출력은 이 표를 쓰지 않음 |
| 검증용 문헌값: Ninković & Adrović (2012, InTech "Gamma Radiation" 1장), Radiological Health Handbook (1970, 미국 PHS), NCRP Report 49 (1976), NUREG/CR-5740 Table 3 | — | `validation.py` 의 비교용 숫자 몇 개(인용) | 각 문헌 저작권 — 숫자 인용·출처 표기 | 각 문헌 | 가능(인용) | 계산에는 쓰지 않고 selftest 비교에만 씀 |
| 원자력안전법 시행령 별표 1 선량한도 | — | `shield.py` `LIMITS` | 법령(저작권 보호 대상 아님) | https://www.law.go.kr | 가능 | 개정 여부는 보건물리 담당자 확인 |
| Google Gemma 4 `gemma4:31b` (Ollama GGUF 빌드) | — | 외부 서비스(OpenAI 호환 API, 자연어 질문 → 입력값·설명) | Apache-2.0 (Ollama `api/show` 의 license 필드로 확인) | https://ai.google.dev/gemma | 가능 | 가중치 동봉 안 함 |
| Ollama | — | 외부 서비스(LLM 서버) | MIT | https://github.com/ollama/ollama | 가능 | 동봉 안 함 |

## 확인 필요 (상업 이용 전)

- **ICRP-74 환산계수 표** ([#1](https://github.com/gggg8657/shield-local/issues/1)) — ICRP 의 표 인용 조건을 확인하지 못했습니다. 다른 공개 프로젝트는 ICRP 유래 계수를 비상업용으로 취급합니다. 회피: 상업 이용 전 ICRP 에 확인하거나, 표를 저장소 밖 파일로 빼고 없으면 공기커마만 표시(제안).
- **축적인자 계수의 전사 경로** ([#2](https://github.com/gggg8657/shield-local/issues/2)) — 값은 공공 문서(NUREG/CR-5740)에서 왔지만 출발점이 비상업 라이선스 프로젝트(Radiographer)의 전사본입니다. 회피: 남은 3 재질(Cu·Al·air)을 스캔 원본과 직접 대조해 독립 전사로 만들기.
- **NIST SRD 8·126 / DDEP 데이터 재배포** ([#3](https://github.com/gggg8657/shield-local/issues/3)) — NIST SRD 는 저작권을 주장하는 데이터이고 DDEP 조건은 미확인입니다. 회피: 데이터 파일을 저장소에서 빼고 설치 때 `tools/build_data.py` 로 받게 하기(제안, 미적용) 또는 NIST 에 확인.
- **파동이(KAERI 캐릭터) — 옛 스크린샷 이력** ([#4](https://github.com/gggg8657/shield-local/issues/4)) — 커밋 715e836 이전의 `docs/img/result.png` 에 포털 펫 파동이가 찍혀 있었습니다. 지금 파일에는 없지만 공개 이력에 남아 있습니다.

이 문서는 법률 자문이 아니며, 라이선스는 각 원문이 우선합니다.

# HANDOFF — 서버에서 바로 실행하기

서버 접속되면 **아래 3줄이면 끝**입니다.

```bash
git clone <이 저장소 주소> && cd Findaptamer
bash scripts/setup_server.sh                                  # 설치 + 수용체 + pool
python src/oracle/calibrate.py --receptor data/mmp9/mmp9_receptor.fasta   # 오라클 검증(게이트)
bash scripts/run_pilot.sh                                     # 설계 실행
```

결과: `results/pilot/top_candidates.csv` (랭킹된 후보 압타머 서열).

---

## 확정된 연구 방향

- **목표: 고친화 결합체(binder)** — 진단/검출용. 억제제가 아니므로 FnII 활성화
  이슈는 해당 없음. 검증은 MST/SPR Kd + MMP2/MMP7 선택성.
- **이길 기준: F3B ≈ 20 nM** (2′-F RNA, MMP9 선택적). 목표는 이에 필적하는 **DNA** 결합체.
- **방법: 오라클 인루프** — AiDTA의 대리보상(2차구조 자기일치도, 단백질 없음)을
  **Boltz-2 공동접힘 결합신호**로 교체. 상세: `docs/02_improved_method.md`,
  `docs/04_boltz2_oracle_ko.md`. AiDTA 비판: `docs/01_aidta_analysis.md`.
- **하드웨어 근거**: `docs/03_hardware_ko.md` (A100 40GB 권장, 24GB가 실용 최소).

## 구현 완료 — 전부 CPU에서 실행 검증됨

| 파일 | 역할 | 검증 |
|---|---|---|
| `src/fragments/build_library.py` | ss 5,440 + ds 5,440 fragment 라이브러리 | ✅ 논문과 동일 |
| `src/scoring/secondary_structure.py` | DNA 접힘 + AiDTA 보상 + 열역학 지표 | ✅ |
| `src/scoring/benchmark_known.py` | F3B 등 기존 압타머 특성화 | ✅ |
| `src/target/mmp9.py` | MMP9 도메인·에피토프·수용체·벤치마크 | ✅ |
| `src/target/fetch_receptor.py` | UniProt/PDB 다운로드 (서버에서 실행) | 네트워크 필요 |
| `src/generator/assembler.py` | AiDTA 조립 semantics 이식 | ✅ |
| `src/generator/pool.py` | fragment pool (generic / HDOCK 기반) | ✅ 44개 생성 |
| `src/oracle/interface.py` | Oracle ABC + CPU 구조 프록시 | ✅ 퇴화 문제 수정 |
| `src/oracle/boltz2.py` | **Boltz-2 결합 오라클** (GPU) | 배선 검증됨 |
| `src/oracle/calibrate.py` | **오라클 검증 게이트** (F3B vs 랜덤, AUROC) | ✅ |
| `src/pipeline/closed_loop.py` | **닫힌 루프 드라이버** | ✅ end-to-end |
| `scripts/setup_server.sh` | 원커맨드 설치 | ✅ 문법·가드 |
| `scripts/run_pilot.sh` | 원커맨드 실행 (캘리브레이션 미통과 시 거부) | ✅ 가드 동작 |

## 안전장치 (중요)

1. **`Boltz2Oracle`는 점수를 지어내지 않습니다.** boltz가 없으면 명확한 에러로 중단.
2. **`run_pilot.sh`는 캘리브레이션 통과 전에는 실행을 거부합니다.**
   실제로 CPU 프록시 오라클은 AUROC 0.167로 FAIL 판정 → 실행 차단됨을 확인.
   즉 오라클이 결합을 구분 못 하면 설계가 시작되지 않습니다.
3. **CPU 드라이런 가능**: `--oracle proxy`로 기계 장치만 점검(결합 정보 없음 명시).

## 서버 접속 후 순서 (상세)

**0. 환경 확인** — setup 스크립트가 자동 출력하지만, 미리 보려면:
```bash
nvidia-smi; nproc; free -g; df -h ~
```

**1. 설치 + 데이터 준비**
```bash
bash scripts/setup_server.sh
```
설치(ViennaRNA/torch/boltz) → MMP9 수용체 다운로드 → fragment pool 생성 → 자체 점검까지 한 번에.
수용체 도메인 선택은 `--domain catalytic`(기본, 107–443) 또는 `catalytic_nofn`
(FnII 제거 → 토큰 수 감소로 co-folding 비용 절감).

**2. 오라클 캘리브레이션 (반드시 통과해야 함)** ⭐
```bash
python src/oracle/calibrate.py --receptor data/mmp9/mmp9_receptor.fasta
```
알려진 결합체(F3B)가 랜덤/셔플/폴리A보다 높은 점수를 받는지 AUROC로 판정.
**PASS(AUROC ≥ 0.75)** 여야 다음 단계로 갑니다. FAIL이면 수용체 범위(도메인 선택),
에피토프, Boltz 설정을 조정해 재시도.

**3. 설계 실행**
```bash
bash scripts/run_pilot.sh
# 규모 조절: ROUNDS=12 PER_ROUND=128 bash scripts/run_pilot.sh
```

**4. 결과 확인**
```bash
column -s, -t results/pilot/top_candidates.csv | head -20
```

## 남은 단계 (Stage 4–5, 서버에서 추가 구현)

- **합의(consensus)**: 상위 후보를 다른 모델(Chai-1)·다중 시드로 재접힘 → 일치하는 것만 통과.
  ※ AF3는 80GB 필요 + 비상업 라이선스라, 40/48GB 환경에서는 Chai-1 권장.
- **에너지 계층**: 짧은 MD(OpenMM) + MM/GBSA로 ΔG 정량 랭킹.
- **최종**: 10–30개 합성 → MST/SPR로 Kd 측정, MMP2/MMP7 선택성 확인.

## 주의사항

- Boltz 최초 실행 시 **가중치 수 GB 다운로드** — 서버 외부 네트워크 필요.
- `--use_msa_server`는 외부 MSA 서버를 씁니다. 폐쇄망이면 단백질 MSA를 미리 만들어
  `Boltz2Config(precomputed_msa=...)`로 넘기세요(수용체가 고정이라 1회만 만들면 됨).
- 공용 클러스터(Slurm)면 `run_pilot.sh`를 `sbatch` 스크립트로 감싸 제출하세요.

## git 참고

최초 커밋 1개가 개인 이메일로 되어 GitHub에서 "Unverified"로 표시됩니다
(force-push가 샌드박스 정책상 차단). 내용은 정상이며 이후 커밋은 모두 정상입니다.

# HANDOFF

## ⚠️ 먼저 읽을 것 — 현재 상태 (2026-08-01 갱신)

**GPU 파이프라인은 완성됐고, 실행됐고, 그 결과 "이 방법은 안 된다"가 확정됐습니다.**
아래 Slurm 가이드는 여전히 유효하지만, **지금 해야 할 일은 GPU에 있지 않습니다.**

| 실행 | 결과 |
|---|---|
| 오라클 캘리브레이션 (job 3406) | **FAIL** — 신호/잡음 0.03, 20 nM 결합체가 무작위보다 낮음 |
| 포즈 재현성 (job 3407) | **FAIL** — 양성 within-Jaccard 0.0 |
| 정준 대조군 Zif268 (job 3408) | **PASS** — ipTM 0.974±0.0007, 포즈 0.974 |

→ 파이프라인은 정상이고, **공동접힘이 압타머를 못 다룹니다.** (`docs/06`)
→ **캘리브레이션 게이트가 잡음 최적화를 막았습니다.** 이게 설계 의도대로 작동한 것입니다.

**현재 진행 중인 방향** (`docs/07`) — 결합을 예측하지 말고, 이미 실행된 실험적
SELEX에서 **사전확률**을 뽑아 탐색공간을 줄인다. **전부 CPU, GPU 불필요:**

```bash
python src/analysis/g4.py            # 공개 MMP-9 압타머 G4 분석 + TBA 순열검정
python src/design/g4_library.py --n 1000   # 집중 라이브러리 설계
```

핵심 발견 3가지는 `docs/07_g4_discovery_ko.md`에 있습니다. 특히:
**공개된 MMP-9 DNA 압타머는 트롬빈 압타머(TBA)를 품고 있습니다** (p = 0.0035).

---

## 인하대 AIX(Slurm) 클러스터 실행 가이드

이 클러스터는 **Slurm로 자원을 신청해서** 쓰는 구조입니다. 로그인 노드에서 설치를
끝내고(=GPU를 놀리지 않음), GPU는 계산할 때만 잡습니다.

## 실행 3단계

```bash
# 0) 로그인 노드에서: 내 파티션 확인 (백틱 주의)
sacctmgr show assoc format=User,Partition where user=`whoami`

# 1) 로그인 노드에서 설치 (GPU 안 잡음)
git clone <저장소> && cd Findaptamer
bash scripts/setup_login.sh

# 2) GPU 신청 + 실행 (PARTITION만 본인 것으로)
bash scripts/slurm_run.sh          # 파티션 기본값 p2
```
`slurm_run.sh`가 GPU를 잡고 → **오라클 캘리브레이션** → **닫힌 루프 설계**까지 실행합니다.
결과: `results/pilot/top_candidates.csv`

### 옵션
```bash
# 배치(백그라운드) 제출
MODE=batch bash scripts/slurm_run.sh
squeue -u `whoami`;  tail -f logs/mmp9apt-*.out

# 자원/시간 조절 (a6000·a100은 최대 4장, a40은 최대 3장, 시간 최대 7-00:00:00)
GPU_TYPE=a6000 GPU_N=4 TIME=2-00:00:00 bash scripts/slurm_run.sh

# 캘리브레이션만 / 설계만
STAGE=calibrate bash scripts/slurm_run.sh
STAGE=pilot     bash scripts/slurm_run.sh
```

### 시간 제한 대비 (중요)
할당 시간이 끝나도 작업이 사라지지 않습니다. 라운드마다 `checkpoint.json`을 저장하고,
`run_pilot.sh`는 기본으로 **이어서 실행(--resume)** 합니다.
```bash
MAX_HOURS=20 TIME=1-00:00:00 bash scripts/slurm_run.sh  # 20시간 후 안전 종료
# 다음 할당에서 같은 명령을 다시 실행하면 중단 지점부터 이어감
```

## 클러스터 사양 / 정책 요약

| GPU | 메모리 | 1회 최대 | 노드 |
|---|---|---|---|
| **a6000** ⭐ | 48GB | **4장** | sv4ka-n1~n4 |
| a100 | 40GB | 4장 | a100-n1~n4 |
| a40 | 48GB | **3장** | sv8ka-n1~n3 |

- **최대 7일**(`7-00:00:00`). 이 계정(mellab)의 파티션은 **p2**(기본값 설정됨). 메일은 `@inha.ac.kr`만.
- p2 사용 가능 노드: a100-n2~n4, **sv4ka-n2~n4(a6000)**, sv8ka-n2~n3(a40).
- 센터가 사용시간을 **모니터링**하며 정책 위반 시 **작업 취소·패널티**. → 설치는 반드시
  로그인 노드에서, GPU는 계산에만.
- **권장: a6000** — 48GB라 수용체 촉매도메인 337aa를 자르지 않고 그대로 넣을 수 있고,
  4장까지 신청 가능해 처리량이 가장 좋습니다.

## 확정된 연구 방향

- **목표: 고친화 결합체(binder)** — 진단/검출용. 억제제가 아니므로 FnII 활성화
  이슈는 해당 없음. 검증은 MST/SPR Kd + MMP2/MMP7 선택성.
- **이길 기준: F3B ≈ 20 nM** (2′-F RNA, MMP9 선택적). 목표는 이에 필적하는 **DNA** 결합체.
- **방법: 오라클 인루프** — AiDTA의 대리보상(2차구조 자기일치도, 단백질 없음)을
  **Boltz-2 공동접힘 결합신호**로 교체. 상세: `docs/02_improved_method.md`,
  `docs/04_boltz2_oracle_ko.md`. AiDTA 비판: `docs/01_aidta_analysis.md`.
- **하드웨어 근거**: `docs/03_hardware_ko.md`. 이 클러스터에서는 **a6000(48GB)** 권장.

## 구현 완료 — 전부 CPU에서 실행 검증됨

| 파일 | 역할 | 검증 |
|---|---|---|
| `src/fragments/build_library.py` | ss 5,440 + ds 5,440 fragment 라이브러리 | ✅ 논문과 동일 |
| `src/scoring/secondary_structure.py` | DNA 접힘 + AiDTA 보상 + 열역학 지표 | ✅ |
| `src/scoring/benchmark_known.py` | F3B 등 기존 압타머 특성화 | ✅ |
| `src/target/mmp9.py` | MMP9 도메인·에피토프·수용체·벤치마크 | ✅ |
| `src/target/fetch_receptor.py` | 수용체 추출 (서열 저장소에 포함, 오프라인 동작) | ✅ |
| `src/generator/assembler.py` | AiDTA 조립 semantics 이식 | ✅ |
| `src/generator/pool.py` | fragment pool (generic / HDOCK 기반) | ✅ 44개 생성 |
| `src/oracle/interface.py` | Oracle ABC + CPU 구조 프록시 | ✅ 퇴화 문제 수정 |
| `src/oracle/boltz2.py` | **Boltz-2 결합 오라클** (GPU) | 배선 검증됨 |
| `src/oracle/calibrate.py` | **오라클 검증 게이트** (F3B vs 랜덤, AUROC) | ✅ |
| `src/pipeline/closed_loop.py` | **닫힌 루프 드라이버** (체크포인트/재개) | ✅ end-to-end |
| `scripts/setup_login.sh` | 로그인 노드 설치 (GPU 미사용) | ✅ 문법 |
| `scripts/slurm_run.sh` | GPU 신청 + 실행 (srun/sbatch) | ✅ 문법 |
| `scripts/run_pilot.sh` | 원커맨드 실행 (캘리브레이션 미통과 시 거부) | ✅ 가드 동작 |

## 안전장치 (중요)

1. **`Boltz2Oracle`는 점수를 지어내지 않습니다.** boltz가 없으면 명확한 에러로 중단.
2. **`run_pilot.sh`는 캘리브레이션 통과 전에는 실행을 거부합니다.**
   실제로 CPU 프록시 오라클은 AUROC 0.167로 FAIL 판정 → 실행 차단됨을 확인.
   즉 오라클이 결합을 구분 못 하면 설계가 시작되지 않습니다.
3. **CPU 드라이런 가능**: `--oracle proxy`로 기계 장치만 점검(결합 정보 없음 명시).

## 오라클 캘리브레이션 (반드시 통과해야 함) ⭐

```bash
python src/oracle/calibrate.py --receptor data/mmp9/mmp9_receptor.fasta
```
알려진 결합체 F3B가 랜덤/셔플/폴리A보다 높은 점수를 받는지 **AUROC**로 판정.
**PASS(AUROC ≥ 0.75)** 여야 `run_pilot.sh`가 실행됩니다. FAIL이면 수용체 범위
(`--domain catalytic` 337aa ↔ `catalytic_nofn` 162aa), 에피토프, Boltz 설정을 조정해 재시도.

결과 확인:
```bash
column -s, -t results/pilot/top_candidates.csv | head -20
```

## Stage 4-5 (구현 완료) — 정밀 검증 및 최종 후보

파일럿이 끝나면 이어서 실행합니다.

```bash
bash scripts/run_refine.sh
# 옵션: TOP=20 SEEDS=3 USE_CHAI=1 ENERGY_N=10 NS=1.0 bash scripts/run_refine.sh
```

| 단계 | 파일 | 하는 일 |
|---|---|---|
| **Stage 4 합의** | `src/pipeline/consensus.py` | 후보마다 **다중 시드 Boltz-2**(+선택적 **Chai-1**)로 재예측. 신뢰도 평균/편차와 **접촉 잔기 집합의 Jaccard 일치도**를 계산 → *"매번 같은 자리에 붙는가"* 를 검증. 단발성 고득점으로는 통과 불가 |
| **Stage 5 에너지** | `src/pipeline/energy.py` | 짧은 **MD(OpenMM, 암시적 용매)** 후 **단일궤적 MM/GBSA**로 ΔG 계산. AiDTA에는 아예 없는 단계 |
| **최종 리포트** | `src/pipeline/report.py` | 전 단계를 합쳐 **합성 후보 shortlist** 생성 (`results/report/REPORT.md`) |

**핵심 지표 읽는 법**
- `contact_jaccard` — 독립 예측 간 접촉 부위 일치도. **낮으면 점수가 높아도 위양성**입니다.
- `dG_bind_kcal_mol` — 엔트로피 무시·암시적 용매이므로 **Kd 예측값이 아니라 순위 신호**입니다.
- `fnII_fraction` — FnII 엑소사이트 접촉 비율. binder 목표에선 참고용이지만, 나중에
  억제제로 전환하면 **역선별 기준**이 됩니다(핵산이 여기 붙으면 MMP9를 활성화).
- `tier` — 3=에너지까지 통과(가장 신뢰), 2=합의까지, 1=파일럿만.

**최종 산출물**
```
results/consensus/consensus_ranked.csv   합의 검증 통과 후보
results/energy/energy_ranked.csv         MM/GBSA ΔG 순위
results/report/REPORT.md                 합성용 shortlist ← 이걸 보세요
```

## 주의사항

- Boltz 최초 실행 시 **가중치 수 GB 다운로드** — 서버 외부 네트워크 필요.
- `--use_msa_server`는 외부 MSA 서버를 씁니다. 폐쇄망이면 단백질 MSA를 미리 만들어
  `Boltz2Config(precomputed_msa=...)`로 넘기세요(수용체가 고정이라 1회만 만들면 됨).
- Boltz 가중치는 `setup_login.sh`가 로그인 노드에서 미리 받습니다(GPU 낭비 방지).
- Stage 5는 **OpenMM 필요**: `pip install openmm`. 미설치 시 자동으로 건너뛰고
  합의 결과만으로 리포트를 만듭니다.
- Stage 4에 Chai-1을 쓰려면 `pip install chai_lab` (Apache-2.0, 상업 가능).

## git 참고

최초 커밋 1개가 개인 이메일로 되어 GitHub에서 "Unverified"로 표시됩니다
(force-push가 샌드박스 정책상 차단). 내용은 정상이며 이후 커밋은 모두 정상입니다.

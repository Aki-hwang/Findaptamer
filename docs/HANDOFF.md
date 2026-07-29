# HANDOFF — 인하대 AIX(Slurm) 클러스터 실행 가이드

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
PARTITION=p1 bash scripts/slurm_run.sh
```
`slurm_run.sh`가 GPU를 잡고 → **오라클 캘리브레이션** → **닫힌 루프 설계**까지 실행합니다.
결과: `results/pilot/top_candidates.csv`

### 옵션
```bash
# 배치(백그라운드) 제출
PARTITION=p1 MODE=batch bash scripts/slurm_run.sh
squeue -u `whoami`;  tail -f logs/mmp9apt-*.out

# 자원/시간 조절 (a6000·a100은 최대 4장, a40은 최대 3장, 시간 최대 7-00:00:00)
PARTITION=p1 GPU_TYPE=a6000 GPU_N=4 TIME=2-00:00:00 bash scripts/slurm_run.sh

# 캘리브레이션만 / 설계만
PARTITION=p1 STAGE=calibrate bash scripts/slurm_run.sh
PARTITION=p1 STAGE=pilot     bash scripts/slurm_run.sh
```

### 시간 제한 대비 (중요)
할당 시간이 끝나도 작업이 사라지지 않습니다. 라운드마다 `checkpoint.json`을 저장하고,
`run_pilot.sh`는 기본으로 **이어서 실행(--resume)** 합니다.
```bash
MAX_HOURS=20 PARTITION=p1 TIME=1-00:00:00 bash scripts/slurm_run.sh  # 20시간 후 안전 종료
# 다음 할당에서 같은 명령을 다시 실행하면 중단 지점부터 이어감
```

## 클러스터 사양 / 정책 요약

| GPU | 메모리 | 1회 최대 | 노드 |
|---|---|---|---|
| **a6000** ⭐ | 48GB | **4장** | sv4ka-n1~n4 |
| a100 | 40GB | 4장 | a100-n1~n4 |
| a40 | 48GB | **3장** | sv8ka-n1~n3 |

- **최대 7일**(`7-00:00:00`), 파티션 `p1/p2/p3`(사용자별 지정), 메일은 `@inha.ac.kr`만.
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

## 남은 단계 (Stage 4–5, 서버에서 추가 구현)

- **합의(consensus)**: 상위 후보를 다른 모델(Chai-1)·다중 시드로 재접힘 → 일치하는 것만 통과.
  ※ AF3는 80GB 필요 + 비상업 라이선스라, 40/48GB 환경에서는 Chai-1 권장.
- **에너지 계층**: 짧은 MD(OpenMM) + MM/GBSA로 ΔG 정량 랭킹.
- **최종**: 10–30개 합성 → MST/SPR로 Kd 측정, MMP2/MMP7 선택성 확인.

## 주의사항

- Boltz 최초 실행 시 **가중치 수 GB 다운로드** — 서버 외부 네트워크 필요.
- `--use_msa_server`는 외부 MSA 서버를 씁니다. 폐쇄망이면 단백질 MSA를 미리 만들어
  `Boltz2Config(precomputed_msa=...)`로 넘기세요(수용체가 고정이라 1회만 만들면 됨).
- Boltz 가중치는 `setup_login.sh`가 로그인 노드에서 미리 받습니다(GPU 낭비 방지).

## git 참고

최초 커밋 1개가 개인 이메일로 되어 GitHub에서 "Unverified"로 표시됩니다
(force-push가 샌드박스 정책상 차단). 내용은 정상이며 이후 커밋은 모두 정상입니다.

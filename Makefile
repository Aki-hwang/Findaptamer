.PHONY: help setup verify run-cpu run-gpu clean
PY ?= python3

help:
	@echo "Findaptamer targets:"
	@echo "  make setup     - create .venv, install CPU stack (+Boltz-2 if a GPU is present)"
	@echo "  make verify    - run the CPU regression checks"
	@echo "  make run-cpu   - closed loop with the CPU StructureProxyOracle -> results/candidates.json"
	@echo "  make run-gpu   - full GPU runbook (Boltz-2 binding oracle) -> results/candidates_boltz2.json"
	@echo "  make clean     - remove venv, __pycache__, generated CSV pools"

setup:
	bash scripts/setup_env.sh

verify:
	$(PY) src/fragments/build_library.py           >/dev/null && echo "[ok] fragment library"
	$(PY) src/scoring/benchmark_known.py            >/dev/null && echo "[ok] benchmark folding"
	$(PY) src/oracle/interface.py                   >/dev/null && echo "[ok] proxy oracle"
	$(PY) src/generator/pool.py                     >/dev/null && echo "[ok] fragment pool"
	$(PY) src/generator/assembler.py                >/dev/null && echo "[ok] assembler"
	$(PY) src/pipeline/closed_loop.py --rounds 4 --pop 24 --out /tmp/findaptamer_verify >/dev/null && echo "[ok] closed loop"

run-cpu:
	$(PY) src/pipeline/closed_loop.py --rounds 12 --pop 60 --top 20

run-gpu:
	bash scripts/run_on_gpu.sh

clean:
	rm -rf .venv
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
	rm -f data/fragments/fragments_ss.csv data/fragments/fragments_ds.csv

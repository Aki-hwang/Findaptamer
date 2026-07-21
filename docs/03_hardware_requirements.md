# Hardware requirements — why the Boltz-2 oracle needs a ≥24 GB GPU

> 한국어 버전: [`docs/03_hardware_requirements.ko.md`](03_hardware_requirements.ko.md)


This note derives, from first principles, why the in-loop binding oracle
(Boltz-2 co-folding MMP9 + a DNA aptamer) does **not** fit on a 4 GB GPU such as
the NVIDIA **T400**, and states the minimum hardware that does. It is meant to be
checkable: every number is either measured, a stated architectural constant, or a
formula you can re-evaluate.

> TL;DR — Co-folding memory is dominated by an **all-pairs (N×N) representation**
> and **triangle operations that scale like N³ in transient compute**. For our
> complex (N ≈ 200–370 tokens) the live activations run to **tens of GB**, and
> the fixed cost (CUDA context + weights + Windows display) already consumes most
> of a 4 GB card before any computation. A 4 GB T400 misses the practical floor
> (~24 GB) by roughly **6–10×**. Even if it fit, its ~1.6 FP32-TFLOPS throughput
> makes a loop of hundreds of co-folds take weeks. Conclusion: a ≥24 GB
> data-center/workstation GPU is required; 4 GB is not a "slow but works" case, it
> is a "does not run" case.

---

## 1. The measured starting point

From the user's workstation (`nvidia-smi`, 2026-07-21):

```
NVIDIA-SMI 595.79   Driver 595.79   CUDA 13.2
GPU 0: NVIDIA T400 4GB   WDDM   657MiB / 4096MiB used   P8
```

Facts this establishes:

| Property | Value | Consequence |
|---|---|---|
| Total VRAM | 4096 MiB (4 GB) | hard ceiling |
| Already in use | 657 MiB | it is also the **display** GPU (WDDM) |
| Driver model | WDDM (not TCC) | extra OS reservation + less predictable free memory |
| Class | Turing TU117, 384 CUDA cores, no usable tensor-core path for this | ~1.6 FP32 TFLOPS |

So the *usable* budget for a computation is **well under 4 GB** — realistically
~2.5–3.4 GB free, and WDDM can claw back more under desktop load.

---

## 2. Why co-folding models are memory-bound (the mechanism)

Boltz-2 and AlphaFold3 are **pairformer / triangle-attention** architectures.
Unlike a sequence model that carries an O(N) hidden state, they carry an explicit
representation of **every pair of tokens**, where a token is one protein residue
or one nucleotide.

Let:

- `N`  = total tokens (protein residues + DNA nucleotides)
- `c_z` = pair-representation channels ≈ **128**
- `c_m` = MSA channels ≈ **256**
- `S`   = number of MSA rows kept (subsampled), often a few hundred
- `H`   = attention heads ≈ **4**
- `b`   = bytes per element (fp32 = 4, bf16/fp16 = 2)

Three structures drive memory:

### 2.1 Pair representation — stored, O(N²)
A tensor of shape `[N, N, c_z]`. One copy:

```
mem_pair = N² · c_z · b
```

### 2.2 Triangle attention — transient, O(N³)
The pair rep is updated by *triangle attention*: for pair (i,j) the update
attends over a third index k. The naive attention-logit tensor is
`[H, N, N, N]`:

```
mem_tri = H · N³ · b
```

This N³ term is the violent one. Real implementations (Boltz included) **chunk**
this so the whole `[N,N,N]` is not materialized at once — trading memory for
time — but even chunked, several large `[N,N,·]` intermediates are live per
block, and there are **tens of pairformer blocks stacked**, each keeping inputs
for the next.

### 2.3 MSA stack + diffusion — additional multipliers
- MSA representation `[S, N, c_m]`, with several live copies in the MSA module.
- The structure is produced by a **diffusion** module that samples multiple
  structures (default several) and runs **recycling** iterations — each multiplies
  the peak by the sample/recycle count.

Peak VRAM is therefore, schematically:

```
Peak ≈ Fixed(context + weights + display)
     + k1·(N²·c_z·b)            ← many live pair-rep copies across depth
     + k2·(N³-ish triangle work) ← chunked, but large and repeated per block
     + k3·(S·N·c_m·b)            ← MSA stack
     + k4·(samples · atoms)      ← diffusion
```

with the `k` multipliers coming from network depth and the number of live
buffers. The key point is **super-linear growth in N**: doubling the complex size
more than quadruples the dominant terms.

---

## 3. Plugging in *our* complex

MMP9 catalytic domain (UniProt P14780, residues 107–443) = **337 aa**; the
FnII-excised "mini" construct = **162 aa**; DNA aptamer ≈ **35 nt**.

| Construct | N (tokens) |
|---|---|
| full catalytic + aptamer | 337 + 35 = **372** |
| FnII-dropped mini + aptamer | 162 + 35 = **197** |

### 3.1 One pair-rep copy `N²·c_z·b`

| N | fp32 | bf16 |
|---|---|---|
| 197 | 19.9 MB | 9.9 MB |
| 372 | 70.8 MB | 35.4 MB |

Small on its own — but there are **many** such copies live across the stack.

### 3.2 One naive triangle-logit tensor `H·N³·b` (fp32, H=4)

| N | H·N³·b |
|---|---|
| 197 | **0.12 GB** |
| 372 | **0.82 GB** |

That is **one** intermediate, in **one** of the two triangle-attention ops, in
**one** of tens of blocks. This is why the *aggregate* peak — even with chunking
and half precision — lands in the **~16–40 GB** range that AF3/Boltz users report
in practice for complexes of this size. (AF3's own guidance cites tens of GB, up
to 80 GB for large nucleic-acid systems; Boltz-2's practical floor for a
protein–DNA co-fold is ~24 GB — this repo's `docs/02_improved_method.md`
hardware table.)

---

## 4. The 4 GB budget reconciliation

Before a single activation is allocated, a 4 GB T400 must already hold:

| Fixed cost | Estimate |
|---|---|
| CUDA + cuDNN context | 0.5 – 1.0 GB |
| Boltz-2 weights (loaded) | ~1 – 2 GB |
| Windows display (WDDM, measured) | 0.66 GB |
| **Fixed subtotal** | **~2.2 – 3.7 GB** |

That leaves **≈ 0.3 – 1.5 GB** for activations, against a requirement of **tens
of GB**. The gap is not a few percent to be tuned away with a flag — it is
**one to two orders of magnitude**.

```
Required (activations)     ~16,000–40,000 MB
Available on T400 after     ~300–1,500 MB
fixed costs
                            ───────────────────
Shortfall                   ~10×–100×
```

Result: `torch.cuda.OutOfMemoryError: CUDA out of memory` during the first
pairformer/triangle pass, for either construct, under default settings.

---

## 5. The second wall: throughput

Suppose, hypothetically, low-memory chunking squeezed the *mini* construct into
memory. Compute then becomes the blocker:

| GPU | FP32 TFLOPS (approx) | Relative |
|---|---|---|
| T400 (Turing TU117) | ~1.6 | 1× |
| RTX 4090 | ~82 | ~50× |
| A100 80 GB | ~19.5 FP32 / ~312 TF32 | ~12×–200× |

A single co-fold that takes ~1–3 minutes on an A100 (MSA cached) would take
**tens of minutes to hours** on a T400 — and the closed loop
(`src/pipeline/closed_loop.py`) issues **hundreds to thousands** of oracle calls
(`rounds × pop`). At, say, 720 evaluations that is **weeks of wall-clock**, with
OOM risk on every call. Not viable.

---

## 6. What a 4 GB card *can* co-fold

To be fair and precise: 4 GB is not useless for these models — it is just far
below *this* target. The N³/N² scaling means only **small** systems fit:

| System size | 4 GB feasibility |
|---|---|
| N ≲ 80–100 (short monomer / small peptide) | often OK with low-mem settings |
| N ≈ 197 (mini MMP9 + aptamer) | unlikely; extreme low-mem only, glacial |
| N ≈ 372 (full MMP9 + aptamer) | no |

Our target sits well past the point where 4 GB stops being an option.

---

## 7. Minimum viable hardware

From `docs/02_improved_method.md` (hardware table), matched to the mechanism above:

| Stage | Minimum | Comfortable |
|---|---|---|
| **Boltz-2 in-loop oracle (Stage 3)** | 1× GPU **24 GB** (RTX 4090 / L40 / A5000) | 1× A100 40 GB |
| AF3 consensus for nucleic acids (Stage 4) | 1× A100/H100 **80 GB** | H100 80 GB |
| MD + MM/GBSA (Stage 5) | RTX 4090 24 GB | A100 |
| System RAM / disk | 64 GB / ~50 GB free (weights+MSA cache) | 128 GB / 1 TB |

Practical ways to obtain a ≥24 GB GPU:

- **Cloud rental** (fastest): RunPod, Vast.ai, Lambda, or Colab (A100). Spin up a
  24–80 GB Linux box; `scripts/run_on_gpu.sh` runs the whole Stage-3 loop.
- **A real workstation GPU**: RTX 4090 / A5000 / A6000, etc.
- On Windows, use **WSL2** (the repo's `.sh`/Makefile need a POSIX shell) or a
  native Linux host. CUDA-on-WSL2 inherits the Windows NVIDIA driver.

---

## 8. How to verify this empirically (don't take it on faith)

On any candidate GPU, before committing to a full run:

```bash
# 1. See the card and its memory
nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv

# 2. Attempt one co-fold and watch memory climb until it either finishes or OOMs
bash scripts/run_on_gpu.sh        # step 3 is a single smoke-test co-fold
# in a second terminal:
nvidia-smi dmon -s m              # live memory monitor
```

On a 4 GB T400 the expected outcome is a `CUDA out of memory` error during the
first triangle-attention pass — the definitive, zero-cost confirmation of this
analysis.

---

## 9. What this does *not* block

The CPU-validated half of the pipeline is unaffected and already committed:
fragment library, DNA folding/scoring, the faithful AiDTA assembler, the
`StructureProxyOracle` closed loop, and the ranked structural-prior shortlist in
`results/candidates.json`. Those run on any machine. Only the **binding** oracle
(Stage 3+) needs the ≥24 GB GPU; it is a one-flag swap (`--oracle boltz2`) once
that hardware is available.

---

### References
- AlphaFold3, *Nature* 2024 (s41586-024-07487-w) — triangle attention; nucleic-acid
  memory guidance.
- Boltz-2, bioRxiv 2025.06.14.659707 (MIT) — co-folding; practical VRAM footprint.
- NVIDIA T400 product spec (Turing TU117, 4 GB GDDR6, 384 CUDA cores).
- This repo: `docs/02_improved_method.md` §"Hardware requirements".

# Frozen CSTR design validation

This audit uses `Lab_Scale_CSTR_Simulation_2_corrected.ipynb`, not the simplified companion simulator. It freezes the 14 schedules in `output/original_schedule_ablation/design_manifest.csv` and reconstructs the same fixed C1 benchmark (f = 0.1) from the near-neutral reference tuning.

## Reproduction

Use Python with NumPy, SciPy, pandas and Matplotlib installed. From the project directory:

```sh
python robustness/run_validation.py
python robustness/refine_validation.py
python robustness/check_implementation.py
python robustness/summarize.py
```

The saved Markdown report is included under `output/robustness/`. Word-document maintenance tools are outside this notebook repository.

## Experimental scope

- 585 numerical design–scenario runs, including exact interval propagation, default/tight RK45, control intervals down to 0.625 s and fractional-delay checks.
- 1,170 unseen runs at 5 s plus the same 1,170 conditions at 2.5 s.
- 18 trajectory comparisons on three selected designs, with implementation parity and reproducibility checks.
- Three new scenarios, 15 frozen designs and 26 nominal/uncertainty conditions. No controller tuning occurs during evaluation.

`protocol.json` was written before the first evaluation; `refinement_protocol.json` records the subsequent numerical follow-up. Protocols record the original simulator and manifest hashes. Raw and paired results are under `output/robustness/`; `summary.txt` contains the compact numerical evidence, and `validation.txt` records implementation checks.

## Interpretation

The primary unseen comparison was specified at 5 s. The 2.5 s results are sensitivity checks, not an opportunity to select a more favourable result. Sampling changes the actual discrete controller and integer-delay representation; it is not solely an ODE integration parameter. Gaussian measurement noise is sampled independently each controller interval. Common seeds produce common draws across designs at the same interval, but the 5 s and 2.5 s noisy cases are not identical continuous-time noise realizations.

Whole-trace IAE includes constant-setpoint disturbances. Recovery intervals are bounded by every setpoint or disturbance event and include the initial interval. An interval that ends outside the ±0.05 pH band has infinite settling time. No failed interval is silently dropped. The test matrix is illustrative and does not define a probability distribution over real plant conditions.

The exact integrator solves the same strong-ion ODE over each held-input interval, retaining the original controller body and event semantics. Fractional delay splits integration at delayed command changes. The stable pH root is algebraically equivalent to the original root. The original simulator source remains unchanged. FOPDT identification is recomputed only to recover the original averaged delay; schedules and tuning selections are loaded from the saved manifest.

The Word report's historical 450 ks wide scenario differs from the corrected notebook's 150 ks development ladder. Historical plots and mean per-step IAE are explicitly separated from new whole-trace unseen results.

Numerical convergence and robustness across the tested uncertainty set were **not established**. The report reports that outcome rather than presenting finite simulation runs as successful control or a universal ranking.

## Section 10 compensation–tracking link

Run `python robustness/compensation_analysis.py` with the same Python/scipy environment. This reads the frozen Section 10 manifest and existing tracking results, reuses the notebook's analytic process-gain function, and writes `compensation_by_variant.csv`, `compensation_tracking.csv`, `compensation_paired_contrasts.csv`, `compensation_curves.csv`, `compensation_summary.md` and `compensation_validation.json` under `output/original_schedule_ablation/`. It does not rerun tuning or tracking simulations. The notebook's Section 10 calls the same function directly on its current in-memory designs/results.

The primary measure exactly retains Section 6's 200-point pH 5–9 max/min proportional loop-gain definition. A 20,001-point grid checks direction sensitivity. Assertions verify stage coverage, preserved original inputs, fixed free parameters and Ti/Td curves, retained Kc anchors, scalar-gain invariance and unchanged proportional compensation under Ti/Td-only resampling. The finding is generated conditionally from isolated anchor pairs, not hard-coded as a successful outcome.

## Central result: loop-gain flatness across tuning methods

Run `python robustness/loop_gain_flatness.py` (NumPy, SciPy, pandas and
Matplotlib required). This reproduces the six frozen common-basis controllers;
C3 uses its recorded selected k=5. Section 6 in the notebook additionally evaluates
the earlier conventional, Rivera, Horn, pure-scheduled and IMC-scheduled designs
using their current in-memory settings. It exports CSV metrics/curves, a Markdown
interpretation and PNG/PDF figures to `output/loop_gain_flatness/`.
The two entry points overwrite that directory with their respective controller sets.

Flatness measures proportional loop gain across pH, with both the legacy 200-point
and refined 20,001-point grids. The normalized curves isolate gain variation from
absolute tuning level. Checks verify scalar invariance, the fixed-controller
baseline and ideal inverse compensation. This central result quantifies static
nonlinearity compensation; Sections 10–11 retain the separate tracking and
robustness evidence. No retuning or new dynamic simulation is performed.

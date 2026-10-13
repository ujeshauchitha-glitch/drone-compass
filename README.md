# drone-compass

Simulation-only feasibility study: **under what conditions, if any, could a radical-pair (spin-chemical) magnetic
compass be useful as an auxiliary heading sensor for GNSS-denied inertial navigation, compared with a conventional
MEMS magnetometer?** A "no, never" answer is a valid outcome.

* No experiment exists. Everything here is simulation. No new physics is claimed.
* Pre-registered in [`PLAN.md`](PLAN.md) before any run. Deviations and bugs are in [`CHANGELOG.md`](CHANGELOG.md).
* Prior art and what is *not* new: [`PRIOR_ART.md`](PRIOR_ART.md). Unverified items: [`open_questions.md`](open_questions.md).

## Layout

| Path | Contents |
|---|---|
| `model/` | spin model (two independent yield implementations), sphere utilities, surrogates, Gate S1 checks |
| `sensor/` | photon-limited K-element sensor, Cramér–Rao bound, Poisson maximum-likelihood axis estimator |
| `nav/` | attitude filters (Stage 3) |
| `experiments/` | one script per stage, writes to `results/` |
| `tests/` | unit tests (`python -m pytest -q`) |
| `results/` | CSV, JSON and figures produced by the experiments |

## Reproduce

```bash
pip install numpy scipy matplotlib reportlab pytest
python -m pytest -q                 # unit tests
python experiments/stage1.py        # Stage 1 (about 8 minutes, one process)
python experiments/stage2_gate.py        # Gate S2 (K3orth, pre-registered) + labelled diagnostic
python experiments/stage2_design.py       # Amendment 1 element-set selection (about 14 min, 2 workers)
python experiments/stage2_sprime.py       # Gate S2', failure rates, Ren comparison (fresh seeds)
python experiments/stage2_noncollinear.py # non-collinear tensor test (about 16 min)
python experiments/stage2_requirements.py # requirements table
python experiments/stage2_failure_threshold.py  # exploratory MC failure-threshold curve (about 14 min)
python experiments/stage2_ren_investigation.py  # CRB-only, seconds
python experiments/stage3_gate.py          # Gate S3 (about 13 min)
python experiments/stage3_grid.py          # Stage 3 main grid (about 36 min, 2 workers)
python experiments/stage3_analysis.py      # decision statistics, break-even, kappa
python experiments/stage3_observability.py; python experiments/stage3_exploratory_snapshot.py
python experiments/stage3_exploratory_tuned.py  # exploratory E1/E2 (about 10 min, 2 workers)
python experiments/stage3_gate_amended.py   # Gate S3(i)' (post hoc, labelled)
python experiments/stage3_cal_grid.py       # calibrated decision study (about 71 min, 2 workers)
python experiments/stage3_cal_analysis.py
python experiments/stage3_tornado.py        # Stage 4 sensitivity (about 27 min, 2 workers)
python experiments/stage3_tornado_analysis.py
python experiments/stage3_grid.py kappa_EXPLORATORY; python experiments/stage3_kappa_analysis.py
python experiments/make_abstract.py; python experiments/make_report.py   # abstract_draft.md, report.pdf
python experiments/make_gate_table.py; python experiments/make_stage2_figures.py; python experiments/make_interim_abstract.py
```

All randomness is seeded from `SeedSequence(20261013)`; runs reproduce exactly. Python 3.14, numpy 2.4, scipy 1.17.

## Status

Stages 0-4 done. **Gates:** S1 passed; S2 (K3orth) FAILED; S2' (amended) passed; S3(i) FAILED, S3(ii) and S3(iii) passed; S3(i)' passed (post hoc, labelled). See `results/gate_table.md`.

**Answer to the pre-registered question: NO** (0 of 390 configurations reduce yaw RMSE by at least 20% at both 60 s and 600 s relative to an uncalibrated MEMS magnetometer; that filter was inconsistent).
**With online calibration (Amendment 3, specified after the NO): 12 of 270** configurations meet the criterion, 11 of them where the MEMS magnetometer is noisy (2 uT) or the dip is 15 deg and the radical-pair sensor reaches 0.3 deg per axis; none at the base setting (0.5 uT, 45 deg).

Deliverables: `report.pdf`, `abstract_draft.md`, `paper_outline.md`, `REVIEW.md`, `PLAN.md` (with Amendments 1-3), `CHANGELOG.md`, `PRIOR_ART.md`, `open_questions.md`, `results/`.

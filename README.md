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
python experiments/make_gate_table.py; python experiments/make_stage2_figures.py; python experiments/make_interim_abstract.py
```

All randomness is seeded from `SeedSequence(20261013)`; runs reproduce exactly. Python 3.14, numpy 2.4, scipy 1.17.

## Status

Stages 0-2 done. Gate S1 passed; Gate S2 (K3orth, pre-registered) FAILED; Gate S2' (amended, selected set, fresh data) passed.
See `results/gate_table.md` and `CHANGELOG.md`. **Stage 3 (navigation) has not been run.**

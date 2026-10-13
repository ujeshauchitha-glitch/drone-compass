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
python experiments/stage2_gate.py   # Stage 2 gate
```

All randomness is seeded from `SeedSequence(20261013)`; runs reproduce exactly. Python 3.14, numpy 2.4, scipy 1.17.

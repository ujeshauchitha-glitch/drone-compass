# paper_outline.md — 6-page PLANS-style outline

Working title: **When could a radical-pair chemical compass help an inertial attitude filter? A simulation-only requirements study**

Target: IEEE/ION PLANS 2027 (full manuscript due 1 Feb 2027; in-person presentation required). Every number below must be taken from `results/` through the generated claims check; no number in this outline is a result by itself.

| Section | Pages | Content | Figures / tables (source) |
|---|---|---|---|
| Abstract | 0.2 | `abstract_draft.md` (generated; ≈ 300 words) | — |
| 1 Introduction | 0.6 | GNSS-denied heading; magnetometer disturbance problem; radical-pair compass as a candidate auxiliary channel; the question and why "no" is a valid answer; contribution = open requirements framework, no new physics, simulation only | — |
| 2 Background and prior art | 0.6 | Spin-chemistry precision limits (Ren 2021, Procopio & Ritz, Cai, Bezchastnov & Domratcheva), cryptochrome evidence, ChemNav, quantum magnetic navigation (anomaly matching); what is not found (no filter-level evaluation; limited search) | Table: prior-art summary (`PRIOR_ART.md`) |
| 3 Models | 1.2 | 3.1 spin yield model, verified tensors, validation checks and the two independent implementations; 3.2 photon-limited sensor, axis-only observability, Cramér–Rao bound and multi-start ML estimator; 3.3 attitude filter (error-state definition, measurement models, sign-invariant residual), disturbance model (shared by both sensors), online calibration | Fig. 1 yield map (`fig_stage1_yield_map.png`); schematic of the chain |
| 4 Validation | 0.8 | Gate table including the failures: S2 FAILED (mirror ambiguity), S2′ (selection rule, fresh data), S3 FAILED (unmodelled acceleration) and the post-hoc S3(i)′; independent second filter; Ren comparison and the unexplained factor ≈ 12; bugs found and fixed (selection objective blind spot, single-start estimator, κ model) | Table: `results/gate_table.md`; Fig. 2 error vs photons (`fig_stage2_error_vs_N.png`) |
| 5 Results | 1.8 | 5.1 sensor requirements (photons per contrast row, failure thresholds, worst-direction caveat, relaxation margin); 5.2 mirror ambiguity and the non-collinear test; 5.3 the pre-registered navigation answer (NO) and why the baseline was inconsistent; 5.4 calibrated baseline and the two regimes where the criterion is met; 5.5 observability and global sign ambiguity; 5.6 sensitivity (tornado) | Fig. 3 failure rates; Fig. 4 yaw RMSE (pre-registered vs calibrated); Fig. 5 tornado; Table: requirements; Table: passing configurations |
| 6 Limits and threats to validity | 0.5 | No experiment; assumed contrast and error model; attitude only; idealised calibration priors; amendment history and forking paths; multiplicity; what was not modelled | `report.pdf` §6 list |
| 7 Conclusion | 0.3 | Honest bottom line: under shared disturbance the radical-pair channel is not useful against a good calibrated MEMS magnetometer; it could matter only if the magnetometer is poor (≈ 2 µT noise) or the dip is low and the sensor reaches ≈ 0.3°/axis; the requirement is a measured contrast, not photons | — |
| References | 0.3 | Verified references only (see `open_questions.md` for items to drop or fix) | — |

Writing rules: use "simulated", "under shot-noise-limited assumptions", "lower bound"; never "demonstrated", "validated" or "real-time". Keep the failed-gate rows in the paper. State the pre-registered answer first and the calibrated study as a later, separately specified analysis.

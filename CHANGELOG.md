# CHANGELOG.md

Bugs, deviations and clarifications, in order. Thresholds in `PLAN.md` are never changed to rescue a result.

## 2026-10-13
* Initial commit: PLAN.md (pre-registration), PRIOR_ART.md, open_questions.md.
* Pre-registration decision (before any run): the hyperfine tensors recalled in the study brief do not match the primary source (Hiscock et al. 2016, main text). The verified values are primary (set "H16"); the recalled values are the labelled sensitivity set "R". Rationale: hard rule 4 (verified-or-ASSUMED).
* Pre-registration decision: Gate S3(i) NEES is evaluated in the correctly modelled nominal condition (no disturbance), because the disturbance is non-white and cannot be consistent in a filter that treats it as white; NEES under disturbance is reported.
* Stage 1 implementation note: `yield_eig` chunk size is chosen from the Hilbert-space size to bound memory (no effect on results).
* Addition (made before looking at any robustness result, labelled EXPLORATORY, no gate depends on it): robustness model **R3** = R2 + one more ¹H on radical A (a_iso = 0.4 mT, d = 0.05 mT) and one more ¹H on radical B (a_iso = 0.3 mT, d = 0.05 mT), all ASSUMED, to show how quickly the minimal model's anisotropy dilutes. Exploratory relaxation sweep uses a meridian scan (φ = 0, π/2; θ ∈ [0, π/2]) plus 1-D polish rather than a full-sphere search, because each Liouville solve is dense (1296²); labelled as scan-based.

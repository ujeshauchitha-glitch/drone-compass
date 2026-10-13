# PLAN.md — pre-registration

Committed before any experiment is run. Thresholds below are never edited after results are seen.
Bugs and any other change go in `CHANGELOG.md` with what/why. Anything marked **ASSUMED** is not measured or verified; it is swept over the stated range.

## 0. Question (fixed)

"Under what conditions, if any, could a radical-pair (spin-chemical) magnetic compass be useful as an auxiliary heading sensor for GNSS-denied inertial navigation, compared with a conventional MEMS magnetometer?"
A "no, never" answer is a valid outcome. **Simulation only. No experiment exists. No new physics is claimed.**

## 1. Global conventions

* Python only: numpy, scipy, matplotlib (+ reportlab for PDF). QuTiP is not installed; the independent cross-check is a Liouville-space solver (see 3.4).
* Randomness: `numpy.random.SeedSequence(20261013)`; each (stage, condition, trial) gets a child via `spawn_key`; runs reproduce bit-for-bit. Max 2 worker processes.
* Units: hyperfine tensors in mT; angular frequency ω = γ_e·(mT value), γ_e = 1.76085963023×10⁸ rad s⁻¹ mT⁻¹ (CODATA 2018, free electron).
* Molecular frame: z = flavin symmetry axis. Field direction n = (sinθ cosφ, sinθ sinφ, cosθ) in the molecular frame.
* Hamiltonian: H/ħ = γ_e B₀ n·(S_A + S_B) + γ_e Σ_k S_A·A_k·I_k (+ optional terms in 3.5). Hyperfine principal axes coincide with molecular x, y, z (**ASSUMED**).

## 2. Stage 0 (done at commit time) — `PRIOR_ART.md`

Read-and-report. Ren et al. photon requirements recorded: 2.29×10⁸ (O1, default), ≈5.3×10⁷ (O2, default), 3.94×10⁵ (optimised a,b,c, O2) for ε_H < 5°. Their bound is information-theoretic, not Cramér–Rao.

## 3. Stage 1 — spin-dynamics model

### 3.1 Primary model
Radical A with two spin-1 nuclei (N5, N10), radical B with no nuclei. B₀ = 50 µT, singlet-born, Haberkorn recombination with k_S = k_T = k, k = 10⁶ s⁻¹ (primary). No exchange, dipolar or relaxation. Dimension 4·9 = 36.

**Tensors — pre-registered decision (made before any run):**
* **Primary set "H16"** (verified against Hiscock et al., PNAS 2016, main text, "Origin of the Spike" section): N5 = (−0.087, −0.100, 1.757) mT, N10 = (−0.014, −0.024, 0.605) mT.
* **Sensitivity set "R"** (recalled in the study brief; **UNVERIFIED**, does not match H16): N5 = (−0.0989, −0.0987, 1.7569) mT, N10 = (−0.0035, −0.0021, 0.4265) mT. Reported in every Stage 1 table, never used for gates.

### 3.2 Yield
Φ_S(n) = (1/M) Σ_{m,n} |⟨m|P_S|n⟩|² k²/(k² + (E_m − E_n)²) over the eigenbasis of H, M = 9 (nuclear dimension), P_S = ¼ − S_A·S_B.
Signal anisotropy ΔΦ = max_n Φ − min_n Φ over the full sphere. Relative anisotropy = ΔΦ/Φ̄ (Φ̄ = sphere mean).

### 3.3 Conditions
* Primary: H16, k = 10⁶ s⁻¹, B₀ = 50 µT.
* Sweeps (full factorial, both tensor sets): k ∈ {10⁵, 10⁶, 10⁷} s⁻¹ × B₀ ∈ {25, 50, 65} µT.
* Map on the sphere (grid + local refinement of extrema), figure of Φ(θ, φ).

### 3.4 Checks (all must pass — Gate S1)
(a) isotropic hyperfine (replace each tensor by its trace/3) gives direction-independent yield, spread < 1e-9;
(b) H = 0 (all terms) gives yield 1 to < 1e-9;
(c) yields within [0, 1] for all grid points;
(d) Φ(n) = Φ(−n) within 1e-9 on ≥ 500 random directions (if it fails, switch to full-sphere everywhere);
(e) two independent implementations agree within 1e-8 on ≥ 200 random directions: eigenbasis formula (3.2) vs. direct Liouville-space Haberkorn solution Φ = k·Tr[P_S (k − L)⁻¹ ρ₀], L = −i[H,·], ρ₀ = P_S/M, solved by dense linear algebra on vec(ρ);
(f) convergence: sphere-grid anisotropy changes by < 1e-6 absolute between 2000 and 8000 Fibonacci points plus Nelder–Mead extremum polish; spherical-harmonic surrogate (3.6) error < 1e-7 absolute on 1000 held-out directions at the chosen degree.
Also reported (not gates): sign-flip of B₀ leaves Φ unchanged; trace/positivity of ρ in the Liouville solver.

**Gate S1 pass condition:** (a)–(f) pass **and** ΔΦ at the primary condition ≥ 1×10⁻³ (absolute). Failing the anisotropy threshold at the primary condition is a reportable negative result and **stops the study at that point** with a report. Anisotropy is reported for every condition regardless.

### 3.5 Robustness (labelled robustness, not gates)
* **R1 (ASSUMED):** primary model + one ¹H on radical B with tensor a_iso·𝟙 + d·diag(−1,−1,2), a_iso ∈ {0, 0.5, 1.5} mT, d ∈ {0.1, 0.5} mT (6 combinations).
* **R2 (ASSUMED):** R1 (a_iso = 0.5, d = 0.5) + one ¹H on radical A, a_iso = 0.8 mT, d = 0.1 mT.
* **Exploratory (added after the external review, labelled EXPLORATORY in the report):** exchange J ∈ {0, ±0.05, ±0.5} mT; axial inter-radical dipolar coupling of strength D ∈ {0, 0.1, 0.3} mT with its axis along molecular z (**ASSUMED**); isotropic spin relaxation rate r ∈ {0, 10⁵, 10⁶, 10⁷} s⁻¹ (Lindblad, both electrons; solved with the Liouville code). Each applied to the primary H16 model.
* Ordering (orientational) disorder is NOT in the spin model; it enters via the contrast parameter ε (4.1). The review noted that real ensembles (e.g. a randomly oriented frozen glass with photoselection) are not fully oriented.

### 3.6 Surrogate
Real spherical-harmonic expansion with even degrees only (Φ(n) = Φ(−n)), fitted by least squares to Φ on a Fibonacci grid; degree chosen as the smallest even L ∈ {8, 12, 16, 20, 24, 32} reaching the 3.4(f) held-out error. Analytic gradients are used in Stage 2.

## 4. Stage 2 — photon-limited sensor and estimator

### 4.1 Sensor
K elements, the same molecule at K orientations R_j (molecular frame vs. body frame). Element j sees axis direction n_j = R_jᵀ n_body. Mean counts
λ_j = N·(1 + ε·(Φ_j − Φ̄)/Φ̄) + D, counts ~ Poisson(λ_j); N = detected signal photons per element per measurement; Φ̄ = sphere mean of Φ; ε ∈ {1, 0.1} (ε = 1 is a perfectly ordered ensemble with the full singlet-yield modulation transferred to the detected intensity; ε < 1 absorbs partial order and incomplete yield-to-signal transfer — **ASSUMED**). D = unmodulated background counts, known.
Element sets: K = 1 (identity); K = 2 (R₁ = I, R₂ = 90° about x); K = 3 orthogonal (cyclic axis permutations); K = 3 optimised (4.4). Baseline sets are fixed here.

### 4.2 Estimator
Poisson maximum likelihood over the axis (unit vectors modulo sign). 4000-point hemisphere Fibonacci grid, best point refined with L-BFGS-B on the tangent parametrisation. Axis-only: errors are measured as min(∠(û, n), ∠(−û, n)). Polarity is not observable (core finding, not a detail).
300 trials per condition, true directions drawn off-grid, uniform on the hemisphere, shared seeds across conditions. N grid: 10⁵ … 10¹¹, 13 log-spaced points.

### 4.3 Cramér–Rao bound
Analytic Poisson Fisher information over the 2-D tangent space at n: J = Σ_j (∇λ_j)(∇λ_j)ᵀ/λ_j with gradients from the surrogate (finite-difference check < 1e-6 relative). Error-norm lower bound sqrt(tr J⁻¹) (RMS), predicted median of a 2-D Gaussian error σ_med = sqrt(ln 2 · tr J⁻¹). Directions where J is singular (rank < 2) are flagged and reported.

### 4.4 Element-orientation design
Compare K = 1, 2, 3 (baselines) and an optimised K = 3 set. Optimiser: `scipy.optimize.differential_evolution` over 3 rotation vectors (seeded, popsize 20, 100 iterations) followed by Nelder–Mead polish; 10 random restarts of Nelder–Mead for comparison. Objective (E-optimality): maximise min over a 2000-point hemisphere grid of λ_min(J). Report mean, median and **worst-case** sqrt(tr J⁻¹) over directions, not only the mean.

### 4.5 Detector realism (all ASSUMED, labelled sweeps)
* Background fraction β = D/N ∈ {0, 1, 10}.
* Additive Gaussian electronics noise σ_e = ρ·√N, ρ ∈ {0, 1, 10}; likelihood variance λ + σ_e²; approximate FIM Σ (∂λ)²/(λ + σ_e²), exact when electronics noise dominates.
* Common-mode multiplicative gain drift g ~ LogNormal(0, σ_g), σ_g ∈ {0, 10⁻³, 10⁻², 10⁻¹}, g estimated as a nuisance parameter (profile likelihood; K ≥ 2 only; K = 1 is non-identifiable and reported as such).
* Lock-in / bandwidth trade-off: measurement duration T = 1/f_u; N = η·R·T with R the detected photon rate and η = 2/π (square-wave demodulation efficiency, ASSUMED); white electronics noise density ν_e (counts·s⁻¹ᐟ²): σ_e = ν_e·√T. Output: attainable error vs update rate f_u ∈ {1, 10, 100} Hz for given (R, ν_e), and the required R, ν_e at 10 Hz.

### 4.6 Comparison with Ren et al.
Evaluate our sensor model with Ren's phenomenological signal Φ(ζ) = a − b·exp[−(ζ−90°)²/(2c²)] (a = 0.5, b = 0.00127, c = 26°; also the optimised row) for K = 3 orthogonal elements aligned with the symmetry axis choice and compute the N giving a 5° CRB-predicted median error. Compare with their 2.29×10⁸ / 5.3×10⁷ / 3.94×10⁵. State all differences (our bound is Cramér–Rao on a 2-D axis; theirs is an information-theoretic bound on a 1-D head direction with r retinal groups). **Any disagreement beyond one order of magnitude is investigated, never tuned.**

### 4.7 Gate S2 (all must pass)
1. At large N — defined as the three largest N-grid points at which the K = 3 orthogonal, ε = 1, no-realism CRB-predicted median error (median over trial directions) is ≤ 1° — the ratio (Monte-Carlo median error)/(median CRB-predicted σ_med at the trial directions) lies in [0.5, 2].
2. The log–log slope of Monte-Carlo median error vs N (same configuration) over grid points whose CRB-predicted median error lies in [0.1°, 3°] is −0.5 ± 0.1.
3. The comparison to Ren et al. (4.6) is documented, with differences listed.

### 4.8 Requirements table
For each condition (K = 3 orthogonal, K = 3 optimised; ε ∈ {1, 0.1}; each realism setting): detected photons N per measurement, photon rate R at 10 Hz (R = N/(η·T)), contrast ε, and tolerable electronics noise ν_e needed for CRB-predicted *median* error (averaged over directions) and worst-case-direction error of 5° and 1°. Each entry is tagged ASSUMED (depends on ε, η, noise assumptions) or model-derived; none is measured.

## 5. Stage 3 — navigation relevance

### 5.1 Filter
Multiplicative EKF on a unit quaternion + gyro bias (error state: 3 attitude + 3 bias). Runs at 10 Hz on delta-angle increments produced from a 100 Hz truth integration (pre-registered simplification; coning ignored; trajectory is slow). Measurements at 10 Hz: accelerometer as gravity direction (R inflated by a lateral-acceleration variance of (0.3 m/s²)² — **ASSUMED**), plus the arm's magnetic channel. 200 Monte-Carlo runs (batched numpy), 600 s.
* Gyro: constant turn-on bias ~ N(0, bias_inst²) + bias random walk whose std over 600 s equals 10 % of bias_inst (ASSUMED) + white noise with ARW; bias_inst ∈ {1, 10, 100} °/h, ARW ∈ {0.1, 0.3, 1} °/√h (full 3×3 grid, all ASSUMED).
* Accelerometer noise 0.05 m/s² per 10 Hz sample (ASSUMED).
* Initial attitude error: tilt 1° 1-σ, yaw 10° 1-σ (ASSUMED); initial bias error equal to the prior.
* Earth field 50 µT; dip ∈ {15°, 45°, 66°} (ASSUMED sweep), base dip = 45°.
* Trajectory (scripted, deterministic given seed): 0–120 s hover with ±3° tilt oscillation at 0.3 Hz; 120–300 s three slow yaw turns of ±90° at ≤ 10°/s; 300–600 s a figure-eight with yaw following heading, tilt ≤ 10°, |a| ≤ 1 m/s².

### 5.2 Arms
A: gyro + accel. B: + MEMS magnetometer (full vector, noise σ_mag ∈ {0.1, 0.5, 2} µT 1-σ per axis; base 0.5 µT, ASSUMED). C: + RP axis-only measurement. D: B + C.
* RP measurement: unit vector û_meas = normalise(u_true,body + e), e tangent Gaussian with isotropic per-axis σ_rp ∈ {0.3°, 1°, 2°, 5°, 10°} (ASSUMED sweep; each σ_rp is linked to the Stage 2 requirement table: "inside the table" means the table lists a required N ≤ N_cap = 10¹⁰ photons per measurement for ε = 1 or ε = 0.1, N_cap swept {10⁸, 10⁹, 10¹⁰}). Update rate 10 Hz.
* Filter handles the sign ambiguity with a sign-invariant residual: r = û_meas − s·û_pred, s = sign(û_meas·û_pred), tangent-plane covariance.
* **Disturbance applies to both sensors identically.** Total field in body frame b_body = S·Rᵀb_earth + d_hi + c(t)·d_mc, with soft-iron S = I + E (symmetric E, entries ~ N(0, 0.01²), applies equally to both sensors, ASSUMED), hard-iron d_hi constant in body frame, d_mc along a fixed body direction scaled by a throttle profile c(t) ∈ [0, 1]; the vector disturbance RMS σ_d ∈ {0, 2, 10} µT is split equally between d_hi and d_mc. Filters inflate R by σ_d²-per-axis (covariance inflation with known σ_d). **The RP sensor is never credited with immunity.** Any sensor-specific advantage is an explicit separate factor κ (exploratory section only; base κ = 1).
* Initial yaw ambiguity test (EXPLORATORY): initial yaw σ = 180°.

### 5.3 Observability
Formal local observability of attitude+bias error state for arms A–D at representative attitudes (rank/singular values of the stacked measurement-Jacobian / observability matrix, with and without gyro excitation), analytic statement of what the axis measurement adds to yaw and tilt versus the vector measurement, and its dependence on dip angle. Sign ambiguity is a global (not local) effect; it is studied by the 180° initial-yaw exploratory runs.

### 5.4 Decision question (pre-registered)
"Is there any region of the assumed parameter space where arm C or D reduces yaw RMSE at t = 60 s **and** t = 600 s by ≥ 20 % relative to arm B, with the radical-pair parameters inside the Stage 2 requirement table?" Evaluated with paired runs (same noise seeds); the criterion is the point estimate ≥ 20 % **and** a paired bootstrap 95 % CI lower bound > 0. The answer may be NO. The full map is reported either way, plus the break-even σ_rp at which arm C's yaw RMSE equals arm B's.

### 5.5 Experimental design (to bound compute)
Base configuration: dip 45°, σ_mag = 0.5 µT, σ0(yaw) = 10°. Full factorial over gyro (3 × 3), σ_d (3), arms, σ_rp (5) at the base configuration. One-at-a-time variations (dip ∈ {15°, 66°}, σ_mag ∈ {0.1, 2}) at the mid gyro setting (10 °/h, 0.3 °/√h) across σ_d, arms and σ_rp. Runtime estimate printed before any run longer than ~5 minutes.

### 5.6 Gate S3
(i) NEES consistency for every arm: time-averaged NEES (6 dof) over the 200 runs inside the two-sided 95 % χ² interval, evaluated in the **correctly modelled nominal condition** (σ_d = 0, no soft-iron, base config) at t ∈ {60, 300, 600} s; consistency under disturbance is reported but is not a gate (it cannot be, because the disturbance is non-white and the filter treats it as white);
(ii) arm A yaw RMSE agrees with the analytic gyro-drift expectation sqrt(σ_b² t² + ARW² t + σ_bw² t³/3 + σ_0²) within a factor of 2 at t = 60 and 600 s for all nine gyro settings (σ_b = bias_inst, σ_bw = 0.1·bias_inst/√600 s, σ_0 = initial yaw σ);
(iii) an independent non-vectorised rotation-matrix EKF with numerical Jacobians reproduces the batched MEKF yaw error trajectory to < 1e-6 rad on 5 runs for arms A–D.

## 6. Stage 4 — sensitivity and write-up
Tornado chart over every ASSUMED parameter (one-at-a-time between its low/high sweep values) on the decision metric; make-or-break parameters listed. Deliverables: PLAN.md, CHANGELOG.md, PRIOR_ART.md, REVIEW.md (10 strongest hostile objections and whether answered), open_questions.md, repo (`model/`, `sensor/`, `nav/`, `tests/`, `experiments/`, README.md), `results/` (CSV + figures), `report.pdf` (5–8 pages), `abstract_draft.md` (≈250–350 words), `paper_outline.md`, optional `interim_abstract.md` (after Stage 2, INTERIM). Every number in the abstract is copied from `results/` and listed in a claims check. Words "demonstrated", "validated" and "real-time" are not used.

## 7. Stage gates and stop rules
Do not start a stage before the previous gate passes. If a check fails: stop, report, do not work around. Negative results reported plainly; exploratory work in its own section.

---

# AMENDMENT 1 — post-hoc, dated 2026-10-13 07:45 (local), written and committed BEFORE any new run

**This amendment is post-hoc.** It was written after Gate S2 failed for K3orth (see CHANGELOG.md, `results/stage2_gate.json`). The study owner approved option A on 2026-10-13 under the conditions below. **No threshold in this plan is changed.** Everything above this line stays exactly as written; K3orth remains recorded as **FAILED Gate S2**.

## A1.1 Reason
Gate S2 compares a Monte-Carlo estimator error with the *local* Cramér–Rao bound. The local bound cannot describe a multimodal likelihood. The tensors N5 and N10 are diagonal in one shared frame, so the molecule has exact π-rotation symmetries and Φ(±nx, ±ny, ±nz) = Φ(n). An element set built from axis permutations (K2, K3orth) inherits a 4-fold mirror ambiguity, up to ~90° apart. The global maximum-likelihood estimate then picks a mirror image about half the time, so the gate was ill-posed for K3orth. The failure is kept and reported as a finding, not removed.

## A1.2 Element-set selection rule (fixed before running anything new)
* **Candidates:** K = 3 elements, R₁ = I (WLOG; a common rotation of all elements is a rotation of the body frame and does not change the objective on the full hemisphere), R₂ and R₃ parametrised by rotation vectors in [−π, π]³ each (6 parameters).
* **Distinguishability of two axes** n, n′ (hemisphere representatives, axis separation ang(n, n′) = arccos|n·n′|) at fixed N: the symmetrised Poisson KL divergence
  D(n, n′) = ½ [KL(n‖n′) + KL(n′‖n)] = ½ Σ_j (λ_j − λ′_j) ln(λ_j/λ′_j), λ_j = N·s_j(n), s_j = 1 + ε(Φ_j − Φ̄)/Φ̄.
* **Objective:** maximise J(R) = min over design-grid pairs with ang(n, n′) > X of D(n, n′), with **X = 10°**, N = 10⁸, ε = 1, β = ρ = σ_g = 0 (D ∝ N for fixed shape, so the maximiser does not depend on N or on ε in this idealised setting).
* **Design directions:** a fixed 600-point hemisphere Fibonacci grid. **Design seeds:** `SeedSequence(20261013, spawn_key=(3, r))` with r = 0, 1, 2 (three independent restarts).
* **Optimiser:** `scipy.optimize.differential_evolution` (popsize 15, maxiter 60, `polish=False`, seeded) followed by Nelder–Mead polish; the restart with the largest J is the selected set. All three restarts' J values are reported.
* **The earlier generic set (K3generic_EXPLORATORY) stays a labelled diagnostic** unless this rule happens to select it; its J and the J of K3orth are reported for comparison.
* **Exploratory comparators (not gates):** the same rule for K = 2 and K = 4 (R₁ = I plus 1 or 3 free rotations).
* **Confirmation of the selected set's J** is recomputed on an independent 3000-point grid and reported.

## A1.3 Design/evaluation separation
* Design uses only the design grid and design seeds above.
* **Evaluation uses fresh data:** 300 true directions drawn uniformly on the hemisphere with `SeedSequence(20261013, spawn_key=(4, 0))`, and noise streams `spawn_key=(4, 1, i)` for N-grid index i. These are different from the Gate S2 seeds (spawn_key prefix 2) and from the design seeds (prefix 3).

## A1.4 Gate S2′ (amended) — same thresholds as Gate S2, applied to the selected set
Gate S2′ conditions 1–3 are word-for-word the conditions of section 4.7 (ratio of Monte-Carlo median error to median CRB-predicted σ_med within [0.5, 2] at the three largest N-grid points whose CRB-predicted median error is ≤ 1°; log–log slope −0.5 ± 0.1 over grid points whose CRB-predicted median error is in [0.1°, 3°]; documented Ren et al. comparison), with ε = 1 and Poisson noise only. Gate S2′ is reported as a **separate row** next to the failed Gate S2 row for K3orth. S2′ passing does not convert S2 into a pass.
* **Reported with Gate S2′ (not gates):** the ambiguity failure rate, defined as the fraction of trials whose axis error exceeds **10°** (= X), versus N, for K3orth, the selected set, and K3generic, at ε = 1 and ε = 0.1 (fresh seeds, `spawn_key=(5, ...)`).
* **Ren comparison (condition 3):** CRB-based N for a 5° median error with the Ren phenomenological signal (a = 0.5, b = 0.00127, c = 26°, and the optimised row a = 0.9, b = 0.00508, c = 10°) using the selected set, compared with their 2.29×10⁸ / ≈5.3×10⁷ / 3.94×10⁵; all differences listed; any disagreement beyond one order of magnitude is investigated, never tuned.

## A1.5 Does the ambiguity depend on a shared principal frame? (new test)
Hypothesis: the mirror ambiguity needs all hyperfine tensors to be diagonal in one common frame (exact D₂ symmetry). Test (all orientations ASSUMED; the principal values of N5 and N10 are the verified H16 values):
* **NC(α):** N5 as in H16; N10 principal frame rotated by α about the axis (1, 1, 0)/√2, α ∈ {0° (control), 5°, 15°, 30°, 60°}; no other change (B₀ = 50 µT, k = 10⁶ s⁻¹).
* **NC-R2(α):** the R2 model of section 3.5 with N10 rotated by α as above and the radical-B proton tensor rotated by α about (0, 1, 1)/√2, α ∈ {0°, 15°, 30°, 60°}.
* **Quantities reported per model:** (i) symmetry violation = max over 500 random directions and the three sign-flip operations (nx→−nx with ny→−ny; ny→−ny with nz→−nz; nx→−nx with nz→−nz) of |Φ(s∘n) − Φ(n)|; (ii) for K3orth, the minimum over test directions of the mirror-pair divergence D(n, s∘n) at N = 10⁸ (s ranges over the three flips; directions with ang(n, s∘n) < X excluded), and the photon number N_mirror = 10 N/D_min at which the weakest mirror pair reaches D = 10 nats (decisive likelihood ratio e¹⁰); (iii) "ambiguity survives" is declared iff D_min at N = 10⁸ is below 10 nats, i.e. N_mirror > 10⁸. Exact symmetry gives N_mirror = ∞.
* The result is labelled as a robustness/mechanism test; no gate depends on it.

## A1.6 Requirements table (section 4.8), with the two additions requested
For the selected K = 3 set and for K3orth (reference), and for ε ∈ {1, 0.1} plus the two **derived** contrasts ε_R2 and ε_relax (the ratio of the model's relative anisotropy aniso_rel to that of the minimal H16 model, from `results/stage1_robustness.csv`, tagged DERIVED/ASSUMED): detected photons N, photon rate R = N/(η T) at 10 Hz (η = 2/π, ASSUMED), tolerable electronics noise ν_e = ρ√(N/T), for a 5° and a 1° target, for both the median-over-directions CRB-predicted σ_med and the worst-case direction (finite maximum over the 2000-point hemisphere grid; infinite → reported as ∞), for the realism settings of section 4.5. **Additions:** (a) for the baseline setting the *ambiguity failure rate* (error > 10°) from Monte Carlo (300 fresh trials) at the listed N, together with the Monte-Carlo median error at that N as a check of the CRB-derived N; (b) the **relaxation margin**: the anisotropy ΔΦ of the underlying spin model divided by the pre-registered 10⁻³ threshold, for the minimal model, R1, R2, R3 and the relaxation sweep (1e5, 1e6, 1e7 s⁻¹), read from the Stage 1 CSVs. Entries are model-derived or ASSUMED; none is measured.

## A1.7 Scope and stop
* **Stage 3 is not started.** The Stage 3 sections of this plan are unchanged and remain "planned".
* After the amendment and Gate S2′, write `interim_abstract.md` from Stages 0–2 only (≈250–350 words), INTERIM-labelled, every number traced to a file in `results/` via a claims-check table generated by a script (`experiments/make_interim_abstract.py`); Stage 3 described as "planned". Forbidden words: "demonstrated", "validated", "real-time". Then **stop and report**.
* A `results/gate_table.csv` / `.md` records every gate row: S1, S2 (K3orth, FAILED), S2′ (amended).
* Runtime estimate for the new runs: about 20 minutes on one process (three DE restarts × three K values dominate); no more than 2 worker processes.

## AMENDMENT 1.1 — dated 2026-10-13 08:05 (local), before any design run

**Reason (found by a smoke test of the objective, no design run, no evaluation run yet):** the A1.2 objective takes the minimum over *pairs of grid points*. The exact mirror image of a grid point is not a grid point, so on the 600-point grid K3orth scores J = 114.4 nats versus 278.1 for the generic set (N = 10⁸): both look "decisive", and the rule is blind to the exact ambiguity it exists to detect (at N = 10⁸ the grid spacing of ~5° is far coarser than the local resolving power).
**Change (the rule's intent is unchanged; thresholds untouched):** the objective is split into two stages.
* **Stage A (search heuristic):** the grid-pair J of A1.2 is used inside the differential-evolution search only.
* **Stage B (selection and reporting objective, J_ref):** for the 40 design directions n with the lowest grid-pair divergence, take the 3 best-matching grid points n′ with ang(n, n′) > X and minimise D(n, n′) over a *continuous* n′ = cosψ·n + sinψ(cosχ·e₁ + sinχ·e₂), ψ ∈ [X, 90°], χ ∈ [0, 2π), with L-BFGS-B (numerical gradient). J_ref = the minimum of these refined divergences. J_ref → 0 for an exactly ambiguous set and is also small where the Fisher information is singular.
* **Selection:** the restart with the largest **J_ref** is the selected set. All restarts' grid-pair J and J_ref are reported; the confirmation uses the 3000-point grid with the same refinement. K3orth and K3generic_EXPLORATORY are scored with J_ref for comparison.
* Everything else in A1.1–A1.7 stands.

---

# AMENDMENT 2 — Stage 3 specification details, dated 2026-10-13 (local), written and committed BEFORE any Stage 3 run

The study owner approved proceeding with Stage 3 ("proceed with Stage 3"). Section 5 above is unchanged; this amendment only fixes details that section 5 left open, so that nothing is decided after seeing results. **No threshold, gate pass condition or decision question of section 5 is changed**; where a clarification touches a gate, the plain-reading value and the stricter reading are both reported.

## A2.1 Frames, conventions, truth
* Nav frame ENU (x east, y north, z up); gravity g = 9.81 m s⁻² along −z; Earth field b_earth = 50 µT·(0, cos dip, −sin dip) (down in the northern hemisphere). R = body→nav rotation. Specific force f_b = Rᵀ(a_nav + 9.81 ẑ).
* **Error state (nav-frame, left-multiplicative):** R_true = exp([δθ]×)·R̂, δb = b_true − b̂. Filter state: quaternion q̂ (storage) + b̂; covariance 6×6 on (δθ, δb). Yaw error = z-component of the rotation vector of R_true R̂ᵀ (rotation about nav-up), wrapped to (−π, π].
* Propagation (10 Hz, dt = 0.1 s) with measured body delta-angle Δθ_m: R̂ ← R̂·exp([Δθ_m − b̂ dt]×); F = [[I, −R̂ dt],[0, I]]; Q = diag(ARW² dt·I, σ_bw² dt·I). Truth delta-angle = rotation vector of R_kᵀR_{k+10} from a 100 Hz truth; measured Δθ_m = Δθ_true + b_true dt + noise with per-axis std ARW·√dt (ARW converted from °/√h; 1 °/√h = (1/60) °/√s).
* Gyro (ASSUMED, sweeps in 5.1): turn-on bias ~ N(0, bias_inst²) per axis plus a bias random walk with σ_bw·√t std, σ_bw = 0.1·bias_inst/√600 s; initial b̂ = 0, P₀(b) = bias_inst²·I. Initial attitude error δθ₀ ~ N(0, diag(1°, 1°, 10°)²) in the nav frame, P₀(θ) = that diagonal.
* **Accelerometer:** measurement z_a = f_b + white noise with per-axis variance 0.05² + 0.3² (m s⁻²)² (the 0.3 m s⁻² is the pre-registered "lateral-acceleration variance"; it is realised in the truth as additional white vibration noise so that the nominal NEES condition is correctly modelled); predicted f̂ = R̂ᵀ(0,0,9.81); H_θ = R̂ᵀ[(0,0,9.81)]×.
* **Trajectory (deterministic, identical for all runs):** t ∈ [0, 600] s. Roll φ(t) = 3° sin(2π·0.3 t), pitch θ(t) = 3° sin(2π·0.3 t + π/2) up to t = 300 s; amplitudes ramp (smoothstep over 300–310 s) to 5° and frequencies change to 0.2 Hz (roll) and 0.13 Hz (pitch, phase 1 rad) for t ≥ 310 s. Yaw: 0 until 130 s; three turns of ±90° (signs +, −, +) starting at 130, 190, 250 s, each lasting 18 s with yaw rate (90°/18 s)(1 − cos 2πu), peak 10°/s; then from 300 s yaw follows the heading of a figure-eight p(τ) = (A sin ωτ, (A/2) sin 2ωτ), A = 5 m, ω = 2π/90 s⁻¹, τ = t − 300, yaw(t) = 90° + (heading_unwrapped(τ) − heading(0)). Linear acceleration a_nav = p̈ (|a| ≤ 0.06 m s⁻² ≪ 1 m s⁻²); attitude is independent of a_nav (e.g. a thrust-vectored or tethered rig), tilts ≤ 5.4°. Throttle c(t) = clip(0.5 + 0.3 sin(2π·0.05 t) + 0.2|ψ̇|/(10°/s), 0, 1).
* Disturbance and sensors (all ASSUMED): total body-frame field b_body = S Rᵀ b_earth + κ(d_hi + c(t) d_mc), S = I + E (E symmetric, each unique entry N(0, 0.01²); applies to the Earth field, equally to both sensors), d_hi and d_mc with independent per-axis N(0, σ_d²/6), so the nominal vector RMS σ_d is split equally; κ = 1 for the MEMS magnetometer and, in all pre-registered arms, for the radical-pair sensor.
  MEMS magnetometer: z_m = b_body + N(0, σ_mag²) per axis (µT). Radical-pair sensor: û_true = b_body/|b_body|; û_meas = s_k·normalise(û_true + e₁t₁ + e₂t₂), tangent errors t ~ N(0, σ_rp²) per tangent axis (σ_rp in rad), random sign s_k = ±1 each measurement, 10 Hz.
* **Filter measurement models.** Arm B: H_θ = R̂ᵀ[b_n]×, R = (σ_mag² + σ_d²/3)·I. Arm C: with u_n = b_earth/|b_earth|, û_pred = R̂ᵀu_n, s = sign(û_meas·û_pred), residual û_meas − s û_pred, H_θ = s R̂ᵀ[u_n]×, R = (σ_rp² + σ_d²/(3|b_earth|²))·I (rad²). Arm D applies B and C sequentially (their disturbance errors are correlated; the filter treats them as independent — reported, not fixed). Joseph-form covariance update, error-state reset after each update (no reset Jacobian).
* **Common random numbers:** all random inputs are standard-normal draws indexed by run (and time) only, generated once from `SeedSequence(20261013, spawn_key=(10, run))`-derived streams and scaled by the swept parameters, so every arm and every parameter setting sees the same noise (paired runs).

## A2.2 Link from σ_rp to the Stage 2 requirement table ("inside the table")
The Stage 3 sweep value σ_rp (deg, per tangent axis) corresponds to a median-error σ_med = sqrt(2 ln 2)·σ_rp = 1.1774 σ_rp. For each contrast row ε ∈ {1 (ideal), 0.1, ε_R2 = 0.0435 (richer assumed model), ε_relax(1e6) = 0.388, ε_relax(1e7) = 0.0111} (values from `results/stage2_requirements_meta.json`): N_req(σ_rp, ε) = max( N_CRB, N_fail(ε) ), where N_CRB = 10⁸·[σ_med,ref(ε)/(1.1774 σ_rp)]² from the median-direction value of `results/stage2_requirements.csv` and N_fail(ε) = the smallest N with ambiguity failure ≤ 5 % from `results/stage2_failure_threshold.json` (ε = 1: 3.16×10⁴; 0.1: 3.16×10⁶; ε_R2: 3.16×10⁷; ε_relax(1e6): 3.16×10⁵; ε_relax(1e7): 3.16×10⁸). A (σ_rp, ε) pair is "inside the requirement table" at cap N_cap ∈ {10⁸, 10⁹, 10¹⁰} photons per element iff N_req ≤ N_cap. The decision question (5.4) is answered **separately for each ε row and each N_cap**, plus a literal overall answer ("yes if any row/cap satisfies it"); the ideal row ε = 1 is labelled an upper bound on what the sensor could do, not a prediction. Limitations stated in advance: isotropic Gaussian per-axis errors replace the direction-dependent CRB (median-direction value used); up to 5 % of epochs may be ambiguity outliers and are not simulated.

## A2.3 Decision statistics
* Yaw RMSE at t = 60 s and t = 600 s over 200 runs (paired). Relative reduction of arm X vs arm B = (RMSE_B − RMSE_X)/RMSE_B. Paired bootstrap (2000 resamples of run indices, `SeedSequence(20261013, spawn_key=(11,))`): the criterion needs the point estimate ≥ 20 % **and** the 2.5th percentile (two-sided 95 % CI lower bound) > 0, at **both** t = 60 s and t = 600 s.
* Break-even σ_rp: for each (gyro setting, σ_d) the σ_rp at which RMSE_C = RMSE_B, by log-linear interpolation of the RMSE ratio over the σ_rp grid (separately at t = 60 s and 600 s); "none in grid" if C is worse (or better) at every grid point.
* **Exploratory only (labelled):** (a) κ ∈ {0.5, 0.1, 0} on the radical-pair sensor's disturbance response (the sensor-specific advantage factor of 5.2) at σ_d ∈ {2, 10} µT, mid gyro, base config, σ_rp ∈ {0.3°, 1°}, arms B, C, D; (b) initial-yaw ambiguity: Monte-Carlo (20000 draws per case) of a snapshot yaw initialisation (TRIAD-style from the first accelerometer and magnetic measurement; for the axis-only sensor the polarity is chosen so the vertical field component points down) and the fraction of runs with |yaw error| > 90°, for dip ∈ {15°, 45°, 66°}, σ_rp ∈ {1°, 5°, 10°}, σ_d ∈ {0, 2, 10} µT, σ_mag = 0.5 µT. No decision depends on either.

## A2.4 Observability (5.3), specification
Numerical local observability of the 6-state error model for arms A–D along the actual truth attitude: the Gramian W = Σ_k Φ(k,0)ᵀ H_kᵀ R_k⁻¹ H_k Φ(k,0) over (i) the 0–120 s hover window and (ii) the 130–148 s yaw-turn window, dips {15°, 45°, 66°}, σ_mag = 0.5 µT, σ_rp = 1°, σ_d = 0. Reported: singular values of W, numerical rank (relative threshold 10⁻¹²), the null direction(s), and the yaw-error standard deviation implied by W⁻¹ (or "unbounded"). Plus the analytic statement of what the axis measurement adds. Sign ambiguity is a global effect and is treated by A2.3(b).

## A2.5 Gate S3 — details (conditions (i)–(iii) of 5.6 unchanged)
* (i) NEES: nominal condition (σ_d = 0, E = 0, base config: dip 45°, σ_mag = 0.5 µT), mid gyro setting (10 °/h, 0.3 °/√h), arms A, B, and C, D at σ_rp = 2°. Statistic: the ensemble-average NEES at each of t ∈ {60, 300, 600} s, 200 runs, 6 dof, compared with the two-sided χ² interval for 6×200 dof divided by 200. **Because 12 tests are made, the gate verdict uses the family-wise 95 % interval (two-sided α = 0.05/12 per test, Bonferroni); the uncorrected 95 % results are reported next to it.** NEES at the other eight gyro settings is reported, not gated.
* (ii) Arm A yaw RMSE vs sqrt(σ_b²t² + ARW²t + σ_bw²t³/3 + σ₀²) with σ_b = bias_inst, σ_bw = 0.1·bias_inst/√600 s, σ₀ = 10°, within a factor of 2 at t = 60 and 600 s for all nine gyro settings (soft-iron and disturbance irrelevant to arm A).
* (iii) Independent non-vectorised filter: rotation matrices (no quaternions), numerical (central-difference) Jacobians, plain Python loop; yaw-error trajectory difference < 10⁻⁶ rad against the batched quaternion MEKF on 5 runs for each of arms A–D (same random inputs).
* If any gate condition fails: stop and report; no work-around.

## A2.6 Compute design and runtime
Base configuration (dip 45°, σ_mag = 0.5 µT, σ₀ = 10°): full factorial over gyro (3 × 3), σ_d ∈ {0, 2, 10} µT, arms A (9 configs), B (27), C (135), D (135) with σ_rp ∈ {0.3°, 1°, 2°, 5°, 10°}; one-at-a-time at the mid gyro setting for dip ∈ {15°, 66°} and σ_mag ∈ {0.1, 2} µT across σ_d, B, C, D. About 440 configurations; each is one batched run of 200 filters × 6000 steps; at most 2 worker processes; the runtime is estimated and printed before any run over 5 minutes.

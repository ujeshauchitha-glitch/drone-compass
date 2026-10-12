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

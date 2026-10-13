# REVIEW.md — the ten strongest objections a hostile reviewer would raise

Status key: **Answered** (addressed by the study), **Partly** (addressed with a stated limit), **Not answered** (cannot be answered by simulation or was not done).
All numbers are from `results/`; none is a measurement.

1. **"There is no experiment. Every key sensor number (contrast ε, per-axis error σ_rp, photon rate) is made up."**
   **Not answered** (by design). The study is a requirements framework: ε and σ_rp are ASSUMED and swept (`PLAN.md` §4.1, §5.2), the requirement table gives photons needed per assumed contrast (`results/stage2_requirements.csv`), and ε is tied to a spin model through derived contrasts (richer-model and relaxation rows). Whether any real material delivers ε ≥ 0.01 with ordered ensembles is exactly what the study cannot say; the closest experiment (Kerpal 2019) was at 120 K.

2. **"The minimal two-nitrogen spin model grossly overstates the anisotropy."**
   **Partly.** True: ΔΦ = 0.206 for the minimal model versus about 0.01 for the assumed richer models R2/R3, and 1.5×10⁻³ with relaxation at 10⁷ s⁻¹ (`results/stage2_margins.csv`). The study reports these and carries them into the requirements through derived ε rows, but all richer couplings are ASSUMED (no verified Trp-side tensors; `open_questions.md` #8) and exchange/dipolar/relaxation are exploratory scans.

3. **"Your magnetometer baseline is a straw man."**
   **Partly.** It was: in the pre-registered main grid the MEMS filter is inconsistent (NEES 183–480 vs ≈ 6) because soft iron is unmodelled, which also produced spurious "wins" at 600 s. Amendment 3 added online hard-iron, soft-iron and motor-current calibration for both sensors; the calibrated filters are consistent (NEES ≈ 6–7) and the baseline improves several-fold. Residual weakness: the calibrated EKF has exact priors, the exact distortion model form and a known throttle profile, which flatters both sensors; a real calibration is harder.

4. **"The pre-registration is not credible: gates failed and you kept amending."**
   **Partly.** Three gates failed (S2 for K3orth, S3(i), and the first-version selection objective was found blind before use). Each failure is kept as a FAILED row (`results/gate_table.md`); every amendment (1, 1.1, 2, 3) was dated and committed before the runs it governs; thresholds were never changed. But S2′ and S3(i)′ were defined after seeing the failure, S3(i)′ in particular after its outcome was visible (labelled post hoc), and the calibrated study was specified after the pre-registered answer (NO) was known. A sceptic may regard the calibrated result as a forking-path result; it is reported next to, not instead of, the pre-registered one.

5. **"Your radical-pair error model is unrealistic: isotropic Gaussian errors, no outliers."**
   **Partly.** The Cramér–Rao bound is direction-dependent; the filter uses a median-direction isotropic value (`PLAN.md` A2.2). The worst field direction needs orders of magnitude more photons (about 4×10⁴× more at ε = 1 for a 5° target, `results/stage2_requirements.csv`), and up to 5 % ambiguity outliers are not simulated. A direction-dependent covariance and an outlier model were not implemented.

6. **"Yaw RMSE of an attitude filter is not GNSS-denied navigation."**
   **Not answered.** Only attitude is estimated; no velocity/position, no drift-vs-time of a position estimate, no visual-inertial coupling. The study argues relevance through heading aiding but does not measure it.

7. **"Both sensors see the same field, so of course disturbance cancels the benefit; a different physical channel might not."**
   **Partly.** The shared-disturbance assumption (κ = 1) is the physically conservative one and no immunity was credited. The sensor-specific factor κ was explored (exploratory, corrected model): the radical-pair channel helps only if it responds to the disturbance at most about half as strongly as a MEMS magnetometer; nothing here supports such a value, and nothing refutes it either.

8. **"Photon counting ignores chemistry: kinetics, bandwidth, temperature, photobleaching, light source."**
   **Not answered.** Only detected photons per measurement at an assumed 10 Hz and lock-in efficiency 2/π are modelled; response time, temperature dependence, ageing and optics are listed under NOT-modelled (`report.pdf` §6). The 120 K experiment of Kerpal et al. is cited as the closest prior art, not matched.

9. **"The scenario is too gentle and the sample too small."**
   **Partly.** The trajectory is slow (tilt ≤ 7°, accelerations ≤ 0.05 m/s², one scripted path), 200 paired runs per configuration with a paired bootstrap, and parameters swept one at a time in the tornado (`results/stage3_tornado_summary.json`); the figure-eight is deliberately benign. No vibration spectrum, accelerometer bias, gyro scale factor, latency or aggressive manoeuvre is modelled, and the unmodelled acceleration already broke one NEES gate.

10. **"270 + 390 configurations, a criterion with a CI, and a 12-of-270 'success' — multiplicity and arbitrary thresholds."**
    **Partly.** The 12 calibrated passes are not isolated: 11 of the 12 cluster in two regimes (MEMS noise 2 µT, or 15° dip) and one D row sits at a high-bias, low-ARW gyro, with effect sizes of 21–60 % and bootstrap lower bounds ≥ 11 %, and the base setting fails. But no multiplicity correction was applied to the 270 tests, the 20 % threshold is a convention, all passes sit at σ_rp = 0.3° except one (the grid edge, so the size of the benefit at smaller σ_rp is unknown), and the bootstrap resamples runs, not scenarios.

**Additional unresolved items (not counted above):** citations that could not be verified (`open_questions.md` #1–#7, e.g. the "Guo et al. 2017" reference was not found), an optimised-signal Ren comparison that remains a factor ≈ 12 off, and a prior-art search limited to web search plus one patent.

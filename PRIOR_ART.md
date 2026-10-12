# PRIOR_ART.md — Stage 0 baseline

Status: written 2026-10-13, before any experiment. Every item is tagged with how it was read:
**[full]** = full text or primary HTML read; **[abstract]** = abstract/metadata only; **[not found]**.
No contents are invented. Anything not opened is listed in `open_questions.md`.

## What is and is not new (stated in every document)

* **Not new:** radical-pair (RP) compass precision limits, Fisher-information-style analyses, optimal and
  robust hyperfine designs, cryptochrome magnetic sensitivity, flight-tested quantum magnetic navigation.
* **New (the only claim):** an open, reproducible, simulation-only requirements-and-feasibility framework that
  carries a spin-chemistry yield model through a photon-limited estimator into an attitude-filter comparison
  against a MEMS magnetometer. No new physics is claimed. **No experiment exists. Everything is simulation.**

## Papers

### Ren, Hiscock & Hore, *Biophys. J.* 120(3):547–555 (2021), "Angular Precision of Radical Pair Compass Magnetoreceptors" [full: PMC7896030 HTML]
* **Shows:** a lower bound on RMS head-direction error ε_H for a retina of r receptor groups, p absorbed photons per group, N = r·p.
  Bound is *information-theoretic* (data-processing inequality → mutual information → Shannon-entropy/variance bound), **not** a Fisher-information / Cramér–Rao bound.
* **Signal model (Eq. 1):** Φ_S(ζ) = a − b·exp[−(ζ−90°)²/(2c²)], ζ = angle between flavin z-axis and the field. Defaults a = 0.5, b = 0.00127, c = 26°.
  Parameters were chosen to approximate a spin simulation of [FAD•⁻ TrpH•⁺] (7 nuclei/radical, 1 µs lifetime, no exchange/dipolar); not fit to experiment.
* **Noise:** per-cell variance σ² = Φ(1−Φ)/p (binomial-like; Poisson not stated explicitly).
* **Photon requirements for ε_H < 5° (Table 2, rp in photons):** default a,b,c, orientation O1: 2.29×10⁸; a=0.9: 8.31×10⁷; +c=10°: 2.63×10⁷; +O2: 6.27×10⁶;
  b=0.00254: 1.57×10⁶; b=0.00508: 3.94×10⁵. (Default with O2 ≈ 5.3×10⁷ from the text.) So "about 3×10⁸" in the study brief ≈ 2.3×10⁸ for O1 and is
  orientation-dependent by ~4×.
* **Axis-only:** model symmetric under ζ → 180°−ζ; Eq. 10 uses |cos ζ|. Polarity is not observable (our reading of the equations; the paper says birds detect the field *axis*).
* **Assumes:** optimal decoder (best case), uniform prior on head direction in horizontal plane, identical cells, axial symmetry, only the [FAD•⁻ TrpH•⁺] pair, 66° inclination geometry.
  A 27-nuclei simulation with exchange/dipolar gave smaller b. Results break down outside 0.1<a<0.9, b>0.1, c<0.1°.
* **Leaves open:** the exact Eq. 5 constants (Hiscock et al. SI needed — not obtained), real devices, ordered ensembles, electronic noise, any navigation use.
* **Consequence for us:** our Cramér–Rao bound is a different quantity from Ren's ε_H; the comparison in Stage 2 is only order-of-magnitude and says so.

### Hiscock et al., *PNAS* 113(17):4634–4639 (2016), "The quantum needle of the avian magnetic compass" [full: PMC4855607 HTML, tensor values read from source text]
* **Verified tensors (principal values, mT; 1 mT ≙ 28 MHz):** N5 = (−0.087, −0.100, 1.757); N10 = (−0.014, −0.024, 0.605). Field 50 µT, singlet yields, singlet-born, first-order rate k with τ = 1/k.
* These **differ from the values recalled in the study brief** (N5 = (−0.0989, −0.0987, 1.7569); N10 = (−0.0035, −0.0021, 0.4265)). N10 differs by 30%. Decision: verified set is primary, recalled set is a labelled sensitivity (see PLAN.md).

### Procopio & Ritz, arXiv:1309.6372, "Robust signals from a quantum-based magnetic compass sensor" [full: arXiv PDF]
* **Shows:** singlet-yield angular sensitivity ΔΦ = Φ_max − Φ_min for one-HF and two-HF models (radical A with k axial hyperfine couplings, radical B none), perfectly ordered crystals, with slow hyperfine fluctuations.
* **Verified numbers:** idealised ΔΦ ≈ 1/8 (region I) and ≈ 1/4 (region II, hyperfine tensor reduced to a vector) in the weak-Zeeman, fast-coherence limit. Their Eq. 2 is the same eigenbasis Lorentzian yield formula used here.
* **Assumes:** perfect order, no dipolar/exchange. **Open:** real molecules, detection, navigation.

### Procopio & Ritz, *J. Chem. Phys.* 152, 065104 (2020), "The reference-probe model for a robust and optimal radical-pair-based magnetic compass sensor" [metadata only; publisher page 403]
* Reference-probe design (one radical without hyperfine coupling). Contents beyond title/authors **not read**; the trade-off between sensitivity and robustness is cited from the arXiv paper above, not from this one.

### Guo et al., *Sci. Rep.* 7, 5826 (2017) [not found]
* Searches did not find a radical-pair compass paper matching this citation. **Not cited as evidence.** Possibly a mis-citation (the Kattnig & Hore 2017 *Sci. Rep.* 7:11640 paper on scavenger amplification is the nearest hit). See `open_questions.md`.

### Xu et al., *Nature* 594, 535 (2021), robin CRY4 [metadata + title]
* In-vitro magnetic sensitivity of European robin CRY4. Biological candidate evidence only; says nothing about device-level heading accuracy. Abstract/full text not read in this study.

### Abrahams et al., "Quantum spin resonance in engineered proteins for multimodal sensing", *Nature* (online 2025; cited as 649, 1172 (2026)) [abstract-level via mirrors]
* Engineered MagLOV proteins show ODMR / radical-pair effects at room temperature in living cells; applications incl. lock-in detection. Whether the work reports direction sensing was **not verified** (brief says none seen in the abstract).

### Cai, *Phys. Rev. Lett.* 106, 100501 (2011), "Quantum probe and design for a chemical compass with magnetic nanostructures" [full: arXiv:1011.5495 PDF text extracted, summary-level]
* Gradient fields (e.g. magnetic nanostructures) can strongly improve directional sensitivity of a chemical compass; suggests a hybrid metallic-organic design. Theory only. Not a navigation evaluation.

### Bezchastnov & Domratcheva, *Sci. Rep.* 16, 1961 (2026) / arXiv:2602.18480 [abstract + intro read; 28-page PDF]
* Arrangement of hyperfine axes and the inter-radical coupling axis can preserve or even sharpen directional sensitivity despite dipolar blurring. The "10–25%" anisotropy figure from the study brief was **not verified** (abstract has no numbers).

### ChemNav, ERC project CORDIS 340451 (Oxford, Dec 2013–Nov 2018) [full: CORDIS final report summary]
* Found: cryptochromes have the needed magnetic-sensor properties; a chemical amplification mechanism (scavenger); maquette proteins with similar field effects; computation of compass precision;
  parallels to polaron pairs in organic semiconductors, with measured magnetoelectroluminescence/magnetoconductance matched to calculated yields.
* **Not found:** any working electronically addressable device or any device-level navigation/heading results. The summary lists no publications; the results tab was not readable.

### Context: quantum magnetic navigation, arXiv:2504.08167 [abstract]
* "Quantum-assured magnetic navigation…" (Muradoglu, … Biercuk): quantum magnetometers + map matching, up to 46× lower positioning error than a strategic-grade INS in airborne trials.
  This is magnetic-*anomaly* navigation with highly sensitive magnetometers, not compassing. The attribution to specific companies in the study brief is **not verified** (author list points to Q-CTRL).

### US12078486B1, Honeywell, "Navigation with magnetic field sensors" [full text read via Google Patents: abstract + claim 1]
* Magnetic gradient sensor + inertial sensors + magnetic field sensor; estimates field change as gradient × inertial velocity and corrects position/velocity. Priority 2021-01-06; granted 2024-09-03.
* No mention of radical pairs, chemical sensing, cryptochrome, or heading from a magnetometer. Overlap with this study: **low** (shared generic idea "magnetic measurement aids an INS error state"). Not a legal opinion.

## Search for navigation-filter-level evaluation of a radical-pair sensor

Searches run 2026-10-13 (web search, standard mode):
1. "radical pair magnetic compass sensor inertial navigation Kalman filter heading aiding chemical compass drone simulation"
2. "cryptochrome-inspired bio-inspired magnetic compass sensor robot navigation extended Kalman filter radical pair"
3. Earlier: reference-probe/Procopio–Ritz, Guo 2017, Cai 2011, ChemNav, US12078486.

**Result: nothing found** that models a radical-pair sensor as an aiding measurement in an IMU/EKF or evaluates it for a drone. This is *absence of evidence from a limited search*, not proof of absence;
patent databases other than the one patent above, IEEE Xplore, Scopus and Google Scholar were not searched.

## What the literature leaves open that this study addresses (simulation only)
1. Photon/contrast/noise requirements for a given axis error, using a spin-chemistry yield model (not a phenomenological one) and a stated estimator.
2. Whether an axis-only (sign-ambiguous) magnetic measurement adds anything to attitude estimation beyond a conventional vector magnetometer.
3. The break-even requirements at which an RP channel would matter in a GNSS-denied attitude filter.

| Gate | Configuration | Status | Detail |
|---|---|---|---|
| S1 | minimal two-nitrogen model, verified H16 tensors | **PASSED** | checks a-f pass; primary anisotropy 0.2057 vs threshold 0.001 |
| S2 (pre-registered) | K3orth, eps=1, Poisson only | **FAILED** | condition 1 ratios 1.448e+04, 2.562e+04, 5.165e+04 (need 0.5-2); condition 2 slope -0.013 (need -0.5 +/- 0.1) |
| S2' (amended, Amendment 1) | K3selected (selected by the pre-registered rule), eps=1, Poisson only, fresh seeds/directions | **PASSED** | condition 1 ratios 0.971, 0.935, 0.937 (need 0.5-2); condition 2 slope -0.523 (need -0.5 +/- 0.1); condition 3: Ren comparison documented in results/stage2_sprime.json |
| S2' diagnostic (not a gate) | K3generic_EXPLORATORY | **DIAGNOSTIC** | ratios 0.928, 0.892, 0.872; slope -0.488 |

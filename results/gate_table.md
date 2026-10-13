| Gate | Configuration | Status | Detail |
|---|---|---|---|
| S1 | minimal two-nitrogen model, verified H16 tensors | **PASSED** | checks a-f pass; primary anisotropy 0.2057 vs threshold 0.001 |
| S2 (pre-registered) | K3orth, eps=1, Poisson only | **FAILED** | condition 1 ratios 1.448e+04, 2.562e+04, 5.165e+04 (need 0.5-2); condition 2 slope -0.013 (need -0.5 +/- 0.1) |
| S2' (amended, Amendment 1) | K3selected (selected by the pre-registered rule), eps=1, Poisson only, fresh seeds/directions | **PASSED** | condition 1 ratios 0.971, 0.935, 0.937 (need 0.5-2); condition 2 slope -0.523 (need -0.5 +/- 0.1); condition 3: Ren comparison documented in results/stage2_sprime.json |
| S2' diagnostic (not a gate) | K3generic_EXPLORATORY | **DIAGNOSTIC** | ratios 0.928, 0.892, 0.872; slope -0.488 |
| S3(i) NEES (pre-registered) | nominal, mid gyro, arms A, B, C and D at sigma_rp=2 deg, t=60/300/600 s | **FAILED** | outside the family-wise 95% interval [5.32, 6.73]: arm A at 600 s: 7.24, arm C at 600 s: 6.83; 10 of 12 tests inside |
| S3(i) diagnostic (not a gate) | same, but with the linear acceleration of the figure-eight removed | **DIAGNOSTIC** | 12 of 12 tests inside the same interval |
| S3(i)' (amended, POST HOC) | same 12 tests on trajectory T0 (linear acceleration = 0); outcome was visible in the diagnostic before the amendment | **PASSED** | 12 of 12 tests inside [5.32, 6.73]; not an independent test |
| S3(ii) arm A vs analytic drift | nine gyro settings, t=60 and 600 s | **PASSED** | ratio simulated/analytic 1.028 to 1.063 (need within a factor 2) |
| S3(iii) independent filter | rotation matrices, numerical Jacobians, 5 runs x arms A-D | **PASSED** | max |yaw difference| 1.19e-08 rad (need < 1e-6) |
| S3 overall (pre-registered) | conditions (i)-(iii) | **FAILED** | (i) failed |
| S3' overall (amended, post hoc) | (i)' and (ii) and (iii) | **PASSED** | reported next to, never instead of, S3 |

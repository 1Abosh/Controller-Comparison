# Loop-gain flatness across tuning methods

F = max(Kc Kp)/min(Kc Kp); F = 1 is perfectly flat. The dense grid contains 20,001 points over pH 5–9. The original 200-point metric is retained separately.

Spread removed = 100 [1 − log(F)/log(F_fixed)]. This compares logarithmic gain variation with the unscheduled plant; negative values mean increased variation. Absolute gain levels are also retained.

| Controller | Dense max/min | Spread (dB) | Spread removed (%) |
|---|---:|---:|---:|
| C1 -- Reference PI fixed | 50.105 | 34.00 | 0.00 |
| C3 -- IMC PI fixed (k=5) | 50.105 | 34.00 | 0.00 |
| C2 -- reference_PI / A_original3 | 72.234 | 37.17 | -9.35 |
| C2b -- reference_PI / E_original5_retuned | 5.046 | 14.06 | 58.65 |
| C4 -- IMC_PID / A_original3 | 72.235 | 37.17 | -9.35 |
| C4b -- IMC_PID / E_original5_retuned | 5.046 | 14.06 | 58.65 |

Central result: quantify how closely each gain schedule compensates the nonlinear local plant gain, while retaining the tuning-dependent absolute loop-gain level. Fixed scalar retuning cannot change flatness; ideal inverse-gain scheduling makes Kc Kp constant. Piecewise interpolation can leave substantial residual variation between anchors.

This is proportional loop-gain flatness across operating pH, not the DC gain of a PI/PID loop (which includes an integrator), nor flatness across frequency. It does not include Ti, Td, sampling, delay, schedule-filter lag or saturation. Use the matched compensation–tracking analysis and unseen robustness results to assess dynamic performance; flatness alone does not establish stability or robust nonlinear control.

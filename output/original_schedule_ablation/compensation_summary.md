# Section 10: compensation paired with tracking

Measure: max(Kc(pH) Kp_local(pH)) / min(Kc(pH) Kp_local(pH)), on the SAME 200-point pH 5–9 grid as Section 6. Lower is flatter. All 14 isolated-family variants are included.

This is static proportional gain compensation, not full PID frequency-response flatness or a stability certificate. Absolute loop-gain levels are retained; a scalar gain change can preserve this ratio while changing tracking.

| Family | Anchor insertion | Ratio before | Ratio after | Scenario | IAE change (%) | Flatter but worse tracking |
|---|---|---:|---:|---|---:|---|
| reference_PI | insert pH 7 only | 21.321 | 7.671 | standardised | +22.17 | True |
| reference_PI | insert pH 8 only | 7.671 | 5.046 | standardised | +12.98 | True |
| reference_PI | insert pH 7 only | 21.321 | 7.671 | wide | +32.04 | True |
| reference_PI | insert pH 8 only | 7.671 | 5.046 | wide | +17.16 | True |
| reference_PI | insert pH 7 only | 21.321 | 7.671 | regulatory | +33.38 | True |
| reference_PI | insert pH 8 only | 7.671 | 5.046 | regulatory | +19.42 | True |
| IMC_PID | insert pH 7 only | 21.321 | 7.671 | standardised | +7.21 | True |
| IMC_PID | insert pH 8 only | 7.671 | 5.046 | standardised | +3.46 | True |
| IMC_PID | insert pH 7 only | 21.321 | 7.671 | wide | +20.15 | True |
| IMC_PID | insert pH 8 only | 7.671 | 5.046 | wide | +19.14 | True |
| IMC_PID | insert pH 7 only | 21.321 | 7.671 | regulatory | +5.25 | True |
| IMC_PID | insert pH 8 only | 7.671 | 5.046 | regulatory | +3.79 | True |

For this benchmark, improved compensation of local process-gain variation did not consistently reduce tracking error at fixed tuning; evaluating schedule refinement therefore requires accounting for its interaction with controller tuning.

12 of 12 isolated anchor/scenario comparisons show both effects. The 20,001-point check preserves every anchor comparison’s compensation direction. Ti/Td curves and existing Kc anchors are verified unchanged in the anchor pairs.

Fixed tuning here means the same family free parameter and entire Ti/Td curves; inserting a Kc anchor deliberately changes the Kc interpolation. These are legacy in-sample benchmark results, not unseen validation. The sampling/convergence limitations in Section 11 still apply. The attribution is conditional on the declared sequence; it is not an order-independent interaction estimate.

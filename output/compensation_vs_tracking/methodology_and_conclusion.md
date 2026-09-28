# Nonlinear compensation versus setpoint tracking

## Methodology

The analysis separates two objectives: compensation of the nonlinear static process gain and dynamic setpoint tracking. It reuses frozen schedules and saved simulation results; no controller is retuned and no new simulation is performed for these figures.

For compensation, evaluate the local proportional loop gain L(pH) = Kc(pH) Kp,local(pH) over pH 5–9. Plot L against operating pH and quantify its flatness using F = max(L)/min(L), where F = 1 denotes perfect compensation. Figure labels use the saved 20,001-point calculation; the plotted curves use the saved 200-point grid. F is one summary over the pH range, not a separate quantity at each pH. It measures proportional gain variation, rather than the full PI/PID frequency response.

For tracking, use IAE = sum(|SP − pH|) Δt, in pH s. Per-step IAE begins at each actual setpoint change and ends immediately before the next change or at the end of the trace. The initial constant-setpoint segment is excluded from per-step metrics. Plot per-step IAE against (a) target pH and (b) the static local loop gain at that target. The latter is interpolated from the saved 200-point curves and is not the time-varying effective gain experienced along a transient. The wide-scenario pH 10 steps are retained in IAE versus pH but omitted from IAE versus local loop gain because the saved gain curves cover only pH 5–9; their exported local gains are blank. S1, S2, etc. identify chronological steps; repeated targets remain separate observations. IAE depends on step direction, size, duration and preceding state, so compare corresponding steps within a scenario rather than inferring a causal pH effect from pooled points.

The primary controlled comparison uses B_analytic3 → C_add7 → D_add8 within each controller family. These stages insert pH 7 and then pH 8 gain anchors while retaining existing Kc anchors, the family free tuning parameter and the complete Ti/Td curves. Kc interpolation deliberately changes. Reference PI and IMC PID are shown separately; standardised and wide setpoint scenarios are also separated. The regulatory scenario addresses disturbance rejection and is excluded from the setpoint-tracking figures. Retuned stages are excluded from this isolated comparison because retuning would confound attribution to anchor insertion.

A fourth figure plots whole-trace IAE against F to connect the global compensation metric directly with dynamic performance. Whole-trace IAE includes the initial interval and any disturbances in the trace; it is distinct from per-step IAE. Lines show the declared schedule-refinement sequence and are not fitted relationships. Absolute local loop gain and flatness are distinct: multiplying all controller gains by a scalar leaves F unchanged but can change IAE.

## Required outputs

1. **Loop-gain flatness versus pH:** `loop_gain_vs_pH.png` shows L(pH), with F in each legend.
2. **IAE versus loop gain:** `IAE_vs_loop_gain.png` shows per-step IAE against local L at the target pH.
3. **IAE versus pH:** `IAE_vs_pH.png` shows per-step IAE against target pH.
4. **Direct conclusion check:** `IAE_vs_flatness.png` shows whole-trace IAE against F.

Each figure also has a PDF version. `per_step_plot_data.csv` preserves all plotted step records. Reproduce with `python robustness/required_outputs.py` in an environment containing NumPy, pandas and Matplotlib.

## Results and conclusion

The isolated anchor sequence reduces the dense-grid flatness ratio from approximately 21.32 to 7.67 and then 5.05 in both controller families. Nevertheless, whole-trace IAE increases at both refinements in both setpoint scenarios: reference PI increases by 22.17% then 12.98% in the standardised scenario and 32.04% then 17.16% in the wide scenario; IMC PID increases by 7.21% then 3.46% and 20.15% then 19.14%, respectively. Thus all eight matched setpoint-scenario comparisons show flatter proportional loop gain alongside higher whole-trace IAE. The separate regulatory comparisons show the same direction but concern disturbance rejection.

For this benchmark and declared refinement sequence, better static nonlinear compensation does not imply better setpoint tracking at otherwise fixed tuning. Schedule evaluation must therefore assess both loop-gain flatness and dynamic IAE, accounting for the interaction with controller tuning. The step-level figures locate the observed tracking errors across transitions; they do not independently establish causation or stability.

These are legacy in-sample results. Existing validation did not establish numerical convergence or robustness across the tested uncertainty set. The figures support a conditional benchmark finding, not a universal claim that flatter gain schedules worsen tracking or that any design is robustly stable.

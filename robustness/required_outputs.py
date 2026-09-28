"""Reproduce compensation versus tracking figures from frozen saved results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'output/original_schedule_ablation'
OUT = ROOT/'output/compensation_vs_tracking'
STAGES = ['B_analytic3', 'C_add7', 'D_add8']
LABELS = ['Analytic 3 anchors', 'Add pH 7', 'Add pH 8']
COLORS = ['#4263eb', '#e67700', '#087f5b']

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    curves = pd.read_csv(SRC/'compensation_curves.csv')
    tracking = pd.read_csv(SRC/'compensation_tracking.csv')
    steps = pd.read_csv(SRC/'original_per_step_metrics.csv')
    families = ['reference_PI', 'IMC_PID']
    scenarios = ['standardised', 'wide']
    steps['target_pH'] = steps.SP_Transition.str.split(' -> ').str[1].astype(float)
    selected = steps[steps.stage.isin(STAGES) & steps.scenario.isin(scenarios)].copy()
    selected['local_loop_gain_at_target'] = [np.interp(r.target_pH,
        curves[(curves.family == r.family) & (curves.stage == r.stage)].pH,
        curves[(curves.family == r.family) & (curves.stage == r.stage)].loop_gain, left=np.nan, right=np.nan)
        for r in selected.itertuples()]
    assert np.isfinite(selected[['IAE', 'target_pH']]).all().all()
    assert selected.loc[selected.target_pH.between(5, 9), 'local_loop_gain_at_target'].notna().all()
    selected.to_csv(OUT/'per_step_plot_data.csv', index=False)
    def save(fig, name):
        fig.tight_layout()
        for ext in ['png', 'pdf']:
            fig.savefig(OUT/f'{name}.{ext}', dpi=180, bbox_inches='tight')
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, family in zip(axes, families):
        for stage, label, color in zip(STAGES, LABELS, COLORS):
            c = curves[(curves.family == family) & (curves.stage == stage)]
            t = tracking[(tracking.family == family) & (tracking.stage == stage)].iloc[0]
            ax.semilogy(c.pH, c.loop_gain, color=color,
                        label=f'{label}; F={t.loop_gain_max_min_dense:.2f}')
        ax.set(title=family, xlabel='Operating pH', ylabel='Proportional loop gain Kc Kp')
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Nonlinear compensation: loop-gain variation across pH')
    save(fig, 'loop_gain_vs_pH')
    for x, name, xlabel in [('target_pH', 'IAE_vs_pH', 'Target pH'),
                            ('local_loop_gain_at_target', 'IAE_vs_loop_gain', 'Local Kc Kp at target pH (interpolated)')]:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for i, family in enumerate(families):
            for j, scenario in enumerate(scenarios):
                ax = axes[i,j]
                for stage, label, color in zip(STAGES, LABELS, COLORS):
                    d = selected[(selected.family == family) & (selected.scenario == scenario) & (selected.stage == stage)]
                    d = d.dropna(subset=[x])
                    ax.scatter(d[x], d.IAE, label=label, color=color, s=35, alpha=.8)
                    for r in d.itertuples():
                        ax.annotate(f'S{r.Step}', (getattr(r,x), r.IAE), xytext=(3,3), textcoords='offset points', fontsize=6)
                if x != 'target_pH': ax.set_xscale('log')
                ax.set(title=f'{family} | {scenario}', xlabel=xlabel, ylabel='Per-step IAE [pH s]')
                ax.grid(alpha=.2); ax.legend(fontsize=8)
        fig.suptitle('Setpoint tracking: matched steps, frozen tuning')
        save(fig, name)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    for ax, family in zip(axes, families):
        for scenario, marker in zip(scenarios, ['o','s']):
            d = tracking[(tracking.family == family) & (tracking.scenario == scenario)].set_index('stage').loc[STAGES]
            ax.plot(d.loop_gain_max_min_dense, d.whole_IAE, marker=marker, label=scenario)
            for stage, r in d.iterrows(): ax.annotate(stage, (r.loop_gain_max_min_dense,r.whole_IAE), fontsize=7)
        ax.set(xlabel='Loop-gain flatness F = max/min (lower is flatter)', ylabel='Whole-trace IAE [pH s]', title=family)
        ax.legend(); ax.grid(alpha=.2)
    save(fig, 'IAE_vs_flatness')
    print(f'Exported four figure pairs and {len(selected)} per-step records to {OUT}')

if __name__ == '__main__': main()

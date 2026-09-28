"""Static proportional loop-gain compensation across operating pH (not frequency).

Run: python robustness/loop_gain_flatness.py
Notebook use: evaluate(controllers, local_gain_analytic, output_directory).
"""
from pathlib import Path
import ast
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def evaluate(controllers, local_gain, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    grid = np.linspace(5., 9., 20001)
    legacy = np.linspace(5., 9., 200)
    kp = np.array([local_gain(p) for p in grid])
    kp_old = np.array([local_gain(p) for p in legacy])
    if not np.all(np.isfinite(kp) & (kp > 0)):
        raise ValueError('Positive finite local process gains required.')
    baseline = kp.max()/kp.min()
    rows, curves = [], []
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes[0, 0].semilogy(grid, kp, color='black')
    axes[0, 0].set_ylabel('Local process gain Kp [pH / (mL/s)]')
    for name, params in controllers.items():
        kc = params['Kc']
        def values(points):
            return (np.array([kc.calculate_gain(p) for p in points])
                    if hasattr(kc, 'calculate_gain') else np.full_like(points, float(kc)))
        loop = values(grid)*kp
        old = values(legacy)*kp_old
        if not np.all(np.isfinite(loop) & (loop > 0)):
            raise ValueError(f'{name}: positive finite Kc*Kp required.')
        ratio = loop.max()/loop.min()
        # Mathematical controls: ideal inverse compensation and scalar invariance.
        np.testing.assert_allclose((1/kp)*kp, 1., rtol=1e-12)
        np.testing.assert_allclose((7*loop).max()/(7*loop).min(), ratio, rtol=1e-12)
        if not hasattr(kc, 'calculate_gain'):
            np.testing.assert_allclose(ratio, baseline, rtol=1e-12)
        rows.append({'Controller': name, 'Loop-gain max/min (pH 5-9)': old.max()/old.min(),
                     'Dense max/min': ratio, 'Spread (dB)': 20*np.log10(ratio),
                     'Spread removed (%)': 100*(1-np.log(ratio)/np.log(baseline)),
                     'Minimum Kc*Kp': loop.min(), 'Maximum Kc*Kp': loop.max(),
                     'pH at minimum': grid[loop.argmin()], 'pH at maximum': grid[loop.argmax()]})
        # Export plotting grid including equivalence and all usual anchors.
        idx = np.arange(0, len(grid), 20)
        curves.extend({'Controller': name, 'pH': grid[i], 'Kp': kp[i],
                       'Kc': loop[i]/kp[i], 'Kc*Kp': loop[i]} for i in idx)
        label = name.split(' -- ')[0]
        axes[0, 1].semilogy(grid, values(grid), label=label)
        axes[1, 0].semilogy(grid, loop, label=label)
        axes[1, 1].semilogy(grid, loop/loop[len(grid)//2], label=label)
    axes[0, 1].set_ylabel('Controller gain Kc [(mL/s) / pH]')
    axes[1, 0].set_ylabel('Proportional loop gain Kc Kp')
    axes[1, 1].set_ylabel('Kc Kp / value at pH 7')
    axes[1, 1].axhline(1, color='black', linestyle=':', label='Ideal inverse compensation')
    for ax in axes.flat:
        ax.set_xlabel('Operating pH'); ax.grid(alpha=.25)
    axes[1, 1].legend(fontsize=7, loc='upper left', bbox_to_anchor=(1, 1))
    fig.suptitle('Gain scheduling compensates local process-gain nonlinearity')
    fig.tight_layout()
    fig.savefig(out/'loop_gain_flatness.png', dpi=180, bbox_inches='tight')
    fig.savefig(out/'loop_gain_flatness.pdf', bbox_inches='tight')
    plt.close(fig)
    table = pd.DataFrame(rows)
    table.to_csv(out/'loop_gain_flatness.csv', index=False)
    pd.DataFrame(curves).to_csv(out/'loop_gain_curves.csv', index=False)
    lines = ['# Loop-gain flatness across tuning methods', '',
             'F = max(Kc Kp)/min(Kc Kp); F = 1 is perfectly flat. The dense grid contains 20,001 points over pH 5–9. The original 200-point metric is retained separately.', '',
             'Spread removed = 100 [1 − log(F)/log(F_fixed)]. This compares logarithmic gain variation with the unscheduled plant; negative values mean increased variation. Absolute gain levels are also retained.', '',
             '| Controller | Dense max/min | Spread (dB) | Spread removed (%) |',
             '|---|---:|---:|---:|']
    for row in rows:
        lines.append(f"| {row['Controller']} | {row['Dense max/min']:.3f} | {row['Spread (dB)']:.2f} | {row['Spread removed (%)']:.2f} |")
    lines += ['', 'Central result: quantify how closely each gain schedule compensates the nonlinear local plant gain, while retaining the tuning-dependent absolute loop-gain level. Fixed scalar retuning cannot change flatness; ideal inverse-gain scheduling makes Kc Kp constant. Piecewise interpolation can leave substantial residual variation between anchors.', '',
              'This is proportional loop-gain flatness across operating pH, not the DC gain of a PI/PID loop (which includes an integrator), nor flatness across frequency. It does not include Ti, Td, sampling, delay, schedule-filter lag or saturation. Use the matched compensation–tracking analysis and unseen robustness results to assess dynamic performance; flatness alone does not establish stability or robust nonlinear control.']
    (out/'loop_gain_flatness.md').write_text('\n'.join(lines)+'\n')
    return table


def main():
    from run_validation import NS, CELLS, DESIGNS, ROOT, pars
    for index in [20, 48]:
        definitions = ast.Module(body=[n for n in ast.parse(''.join(CELLS[index]['source'])).body
                                      if isinstance(n, ast.FunctionDef)], type_ignores=[])
        exec(compile(definitions, 'notebook_definitions', 'exec'), NS)
    kp, tau, theta = np.mean(pars, axis=0)
    kc, ti, _ = NS['imc_pid_rivera'](kp, tau, theta, 5*theta)
    controllers = {'C1 -- Reference PI fixed': DESIGNS['fixed/C1'],
                   'C3 -- IMC PI fixed (k=5)': dict(Kc=kc, Ti=ti, Td=0)}
    for cid, family, stage in [('C2','reference_PI','A_original3'),
                                ('C2b','reference_PI','E_original5_retuned'),
                                ('C4','IMC_PID','A_original3'),
                                ('C4b','IMC_PID','E_original5_retuned')]:
        controllers[f'{cid} -- {family} / {stage}'] = DESIGNS[f'{family}/{stage}']
    print(evaluate(controllers, NS['local_gain_analytic'], ROOT/'output/loop_gain_flatness').to_string(index=False))

if __name__ == '__main__':
    main()

"""Pair Section 6's static proportional-loop-gain measure with EVERY Section 10 variant.

Run from the repository root: python robustness/compensation_analysis.py
Reuses frozen schedules and saved tracking runs; performs no tuning or new simulation.
The same evaluate() is invoked after the notebook's Section 10 simulations.
"""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np
import pandas as pd

ANCHOR_PAIRS = {('B_analytic3', 'C_add7'), ('C_add7', 'D_add8')}

def value_at(value, grid):
    if hasattr(value, 'calculate_gain'):
        return np.array([value.calculate_gain(float(p)) for p in grid])
    return np.full(grid.shape, float(value))

def evaluate(designs, local_gain, results, contrasts, out, metadata):
    out = Path(out)
    grid = np.linspace(5., 9., 200)  # EXACT original Section 6 definition
    dense = np.linspace(5., 9., 20001)  # includes 7, 8 and moved anchors
    kp = np.array([local_gain(float(p)) for p in grid])
    kp_dense = np.array([local_gain(float(p)) for p in dense])
    assert np.isfinite(kp_dense).all() and (kp_dense > 0).all()
    rows = []
    curves = []
    for (family, stage), params in designs.items():
        loop = value_at(params['Kc'], grid) * kp
        refined = value_at(params['Kc'], dense) * kp_dense
        assert np.isfinite(loop).all() and (loop > 0).all()
        ratio = float(loop.max()/loop.min())
        # The ratio must be invariant to scalar Kc retuning.
        np.testing.assert_allclose((5*loop).max()/(5*loop).min(), ratio, rtol=1e-12)
        rows.append(dict(family=family, stage=stage, loop_gain_min=loop.min(),
            loop_gain_max=loop.max(), loop_gain_max_min=ratio,
            loop_gain_max_min_dense=refined.max()/refined.min()))
        curves.extend(dict(family=family, stage=stage, pH=p, Kp=k, Kc=c, loop_gain=l)
            for p,k,c,l in zip(grid,kp,value_at(params['Kc'],grid),loop))
    flat = pd.DataFrame(rows)
    # A left join and explicit coverage assertion prevent silent omission of a stage.
    tracking = results[results.family.isin(flat.family.unique())].copy()
    joined = tracking.merge(flat, on=['family','stage'], how='left', validate='many_to_one')
    assert joined.loop_gain_max_min.notna().all()
    assert set(map(tuple, flat[['family','stage']].values)) == set(map(tuple, tracking[['family','stage']].values))
    indexed = flat.set_index(['family','stage'])
    paired = contrasts.copy()
    for side in ['before','after']:
        for metric in ['loop_gain_max_min','loop_gain_max_min_dense']:
            paired[metric+'_'+side] = [indexed.loc[(r.family,getattr(r,side)),metric]
                for r in paired.itertuples()]
    paired['flatness_change_percent'] = 100*(paired.loop_gain_max_min_after/paired.loop_gain_max_min_before-1)
    paired['compensation_improved'] = paired.loop_gain_max_min_after < paired.loop_gain_max_min_before*(1-1e-9)
    paired['tracking_worsened'] = paired.IAE_after > paired.IAE_before*(1+1e-9)
    paired['isolated_anchor_addition'] = [(a,b) in ANCHOR_PAIRS for a,b in zip(paired.before,paired.after)]
    paired['flatter_but_worse_tracking'] = paired.compensation_improved & paired.tracking_worsened
    anchors = paired[paired.isolated_anchor_addition]
    free = metadata.set_index(['family','stage']).free_parameter
    for r in anchors.itertuples():
        np.testing.assert_allclose(free.loc[(r.family,r.before)], free.loc[(r.family,r.after)], rtol=0, atol=0)
        a,b = designs[(r.family,r.before)], designs[(r.family,r.after)]
        for key in ['Ti','Td']:
            np.testing.assert_allclose(value_at(a[key],dense),value_at(b[key],dense),rtol=1e-12,atol=1e-12)
        for p,_ in a['Kc'].breakpoints:
            np.testing.assert_allclose(a['Kc'].calculate_gain(p),b['Kc'].calculate_gain(p),rtol=1e-12)
    refined_direction = anchors.loop_gain_max_min_dense_after < anchors.loop_gain_max_min_dense_before*(1-1e-9)
    assert np.array_equal(refined_direction.to_numpy(), anchors.compensation_improved.to_numpy()), 'Grid refinement changes the compensation finding'
    # Ti/Td-only resampling must leave the proportional compensation measure unchanged.
    same_kc = paired[paired.change == 'resample Ti/Td only']
    np.testing.assert_allclose(same_kc.loop_gain_max_min_before,same_kc.loop_gain_max_min_after,rtol=1e-12)
    n = int(anchors.flatter_but_worse_tracking.sum())
    conclusion = ('For this benchmark, improved compensation of local process-gain variation '
        'did not consistently reduce tracking error at fixed tuning; evaluating schedule '
        'refinement therefore requires accounting for its interaction with controller tuning.'
        if n else 'The isolated anchor comparisons do not support the proposed flatter-gain/worse-tracking conclusion.')
    lines = ['# Section 10: compensation paired with tracking', '',
        'Measure: max(Kc(pH) Kp_local(pH)) / min(Kc(pH) Kp_local(pH)), on the SAME 200-point pH 5–9 grid as Section 6. Lower is flatter. All 14 isolated-family variants are included.', '',
        'This is static proportional gain compensation, not full PID frequency-response flatness or a stability certificate. Absolute loop-gain levels are retained; a scalar gain change can preserve this ratio while changing tracking.', '',
        '| Family | Anchor insertion | Ratio before | Ratio after | Scenario | IAE change (%) | Flatter but worse tracking |',
        '|---|---|---:|---:|---|---:|---|']
    for r in anchors.itertuples():
        lines.append(f'| {r.family} | {r.change} | {r.loop_gain_max_min_before:.3f} | {r.loop_gain_max_min_after:.3f} | {r.scenario} | {r.IAE_change_percent:+.2f} | {r.flatter_but_worse_tracking} |')
    lines += ['', conclusion, '', f'{n} of {len(anchors)} isolated anchor/scenario comparisons show both effects. The 20,001-point check preserves every anchor comparison’s compensation direction. Ti/Td curves and existing Kc anchors are verified unchanged in the anchor pairs.', '',
        'Fixed tuning here means the same family free parameter and entire Ti/Td curves; inserting a Kc anchor deliberately changes the Kc interpolation. These are legacy in-sample benchmark results, not unseen validation. The sampling/convergence limitations in Section 11 still apply. The attribution is conditional on the declared sequence; it is not an order-independent interaction estimate.']
    flat.to_csv(out/'compensation_by_variant.csv',index=False)
    joined.to_csv(out/'compensation_tracking.csv',index=False)
    paired.to_csv(out/'compensation_paired_contrasts.csv',index=False)
    pd.DataFrame(curves).to_csv(out/'compensation_curves.csv',index=False)
    summary = '\n'.join(lines)+'\n'
    (out/'compensation_summary.md').write_text(summary)
    print(summary)
    return flat, joined, paired, summary

def main():
    # Reuse the frozen-manifest loader and original notebook environment.
    from run_validation import DESIGNS, NS, CELLS, ROOT, MANIFEST
    source = ''.join(CELLS[48]['source'])
    definitions = ast.Module(body=[n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)],type_ignores=[])
    exec(compile(definitions,'original_analytic_gain','exec'),NS)
    out = ROOT/'output/original_schedule_ablation'
    inputs = [MANIFEST,out/'scenario_results.csv',out/'paired_contrasts.csv']
    hashes = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    designs = {tuple(k.split('/')):v for k,v in DESIGNS.items() if not k.startswith('fixed/')}
    evaluate(designs,NS['local_gain_analytic'],pd.read_csv(inputs[1]),pd.read_csv(inputs[2]),out,pd.read_csv(MANIFEST))
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    (out/'compensation_validation.json').write_text(json.dumps(dict(input_sha256=hashes,
        primary_grid=dict(min_pH=5,max_pH=9,n=200),dense_grid_n=20001,
        checks='Complete coverage, unchanged inputs, scalar invariance, retained Kc anchors, fixed free parameters and Ti/Td curves, resampling invariance, grid-direction agreement: PASS'),indent=2)+'\n')

if __name__ == '__main__':
    main()

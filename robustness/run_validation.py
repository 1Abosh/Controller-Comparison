"""Frozen-design validation of the ORIGINAL notebook's incremental PID simulator.

Run with Python + numpy/pandas/scipy/matplotlib. No tuning search is performed.
The original function is loaded and instrumented with checked source substitutions;
its controller, filtering, clipping, event convention and metrics remain unchanged.
"""
from pathlib import Path
import ast
import contextlib
import hashlib
import io
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/robustness'
NOTEBOOK = ROOT / 'Lab_Scale_CSTR_Simulation_2_corrected.ipynb'
MANIFEST = ROOT / 'output/original_schedule_ablation/design_manifest.csv'
CELLS = json.loads(NOTEBOOK.read_text())['cells']
SOURCE = ''.join(CELLS[7]['source'])
NS = {}
with contextlib.redirect_stdout(io.StringIO()):
    for i in [2, 3, 4, 5, 6, 7]:
        exec(''.join(CELLS[i]['source']).replace('import matplotlib.pyplot as plt', ''), NS)
    # Recompute only original identification, not any controller selections.
    pars = [NS['run_step_test'](pH_target=p, step_pct=s, Qb0_guess=g,
                               region_label=str(p), plot=False)
            for p, s, g in [(4, .0125, .05), (6, -.001, 5), (9, .05, 5)]]
NS['theta'] = float(np.mean([p[2] for p in pars]))
NOMINAL = {k: NS[k] for k in ['Qa', 'Ca', 'Cb', 'V', 'Kw']}
FMAX = NS['default_F_max']()
NOMINAL_FLOW = NS['solve_Qb_for_pH']
Scheduler = NS['PiecewiseLinearGainController']

def parse_value(s):
    value = ast.literal_eval(re.sub(r'np\.float64\(([^()]*)\)', r'\1', s))
    return Scheduler(value) if isinstance(value, list) else value

DESIGNS = {}
for row in pd.read_csv(MANIFEST).to_dict('records'):
    DESIGNS[row['family'] + '/' + row['stage']] = {
        k: parse_value(row[k + '_anchors']) for k in ['Kc', 'Ti', 'Td']}
    DESIGNS[row['family'] + '/' + row['stage']]['schedule_filter_tau'] = 60.
# C1 f=.1 is the original near-neutral reference PI; A uses f=.02.
a = DESIGNS['reference_PI/A_original3']
DESIGNS['fixed/C1'] = dict(Kc=5*a['Kc'].calculate_gain(6),
                          Ti=a['Ti'].calculate_gain(6), Td=0.)

def events(times, values, key):
    return [dict(time=t, **{key: v}) for t, v in zip(times, values)]

DEV = {}
exec(''.join(CELLS[56]['source']), NS)
for k in ['standardised', 'wide', 'regulatory']:
    DEV[k] = NS[k + '_scenario']
# Held out from all historical selections; defined before evaluation.
TEST = {
    'reverse_off_anchor': dict(initial_pH=8.7, tend=60000,
        sp_events=events([0, 6000, 24000, 42000], [8.7, 7.35, 5.4, 6.65], 'SP')),
    'cross_equivalence': dict(initial_pH=5.25, tend=60000,
        sp_events=events([0, 7000, 24000, 43000], [5.25, 8.4, 6.8, 7.2], 'SP'),
        dist_events=events([16000, 35000, 51000], [1.03, .97, 1.02], 'Ca_factor')),
    'regulatory_off_anchor': dict(initial_pH=7.35, tend=60000,
        sp_events=events([0], [7.35], 'SP'),
        dist_events=events([6000, 24000, 42000], [1.07, .94, 1.], 'Ca_factor')),
}

def instrument():
    src = SOURCE
    def replace(old, new):
        nonlocal src
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    replace('def run_closed_loop_sim(', 'def audited_sim(')
    replace('title="Closed-Loop pH Control", show_mv=False, plot=True):',
            'title="Closed-Loop pH Control", show_mv=False, plot=True, noise=None, delay_mode="legacy"):' )
    replace('pH_arr[0] = initial_pH',
            'pH_arr[0] = initial_pH\n    measured = np.zeros(n)\n    measured[0] = initial_pH + noise[0]')
    replace('e[0] = sp_at(0) - pH_arr[0]', 'e[0] = sp_at(0) - measured[0]')
    # Replace only feedback references, never the true-state storage or metrics.
    start, end = src.index('    pH_sched_filt ='), src.index('        current_Ca =')
    block = src[start:end].replace('pH_arr[', 'measured[')
    src = src[:start] + block + src[end:]
    replace("        sol = solve_ivp(dxdt_func, (t_span_start, t_span_start + dt), x0,\n                         args=(Qb_delayed, current_Ca), method='RK45')",
        '''        if delay_mode == 'fractional':
            ratio = theta_dead / dt
            lag = int(np.floor(ratio))
            fraction = ratio - lag
            older = F[max(0, k-lag-1)]
            newer = F[max(0, k-lag)]
            mid = advance(dxdt_func, x0, older, current_Ca, fraction*dt)
            sol = advance(dxdt_func, [mid.y[0,-1]], newer, current_Ca, (1-fraction)*dt)
        else:
            sol = advance(dxdt_func, x0, Qb_delayed, current_Ca, dt)''')
    replace('pH_arr[k] = x_to_pH(x0[0])',
            'pH_arr[k] = x_to_pH(x0[0])\n        measured[k] = pH_arr[k] + noise[k]')
    # Stable algebraic evaluation, identical chemistry; audited separately.
    exec(src, NS)
    return src

AUDITED_SOURCE = instrument()
SETTINGS = {}
def advance(fun, x, q, ca, dt):
    if SETTINGS['solver'] == 'exact':
        rate = (NS['Qa']+q)/NS['V']
        ss = (NS['Qa']*ca-q*NS['Cb'])/(NS['Qa']+q)
        value = x[0] + (ss-x[0])*(-np.expm1(-rate*dt))
        return SimpleNamespace(y=np.array([[value]]))
    sol = solve_ivp(fun, (0, dt), x, args=(q,ca), method='RK45',
                    rtol=SETTINGS['rtol'], atol=SETTINGS['atol'])
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol
NS['advance'] = advance
LEGACY_PH = NS['x_to_pH']
def stable_ph(x):
    root = np.sqrt(x*x+4*NS['Kw'])
    h = (x+root)/2 if x >= 0 else 2*NS['Kw']/(root-x)
    return -np.log10(h)

def metrics(res, scenario):
    t, y, sp, u = [res[k] for k in ['t','pH','SP','F']]
    if not all(np.isfinite(v).all() for v in [t,y,sp,u]):
        raise ValueError('Nonfinite trace')
    dt = t[1]-t[0]
    # Include disturbances and initial mismatch, unlike legacy SP-only metrics.
    starts = sorted({0.} | {e['time'] for e in scenario['sp_events']}
                    | {e['time'] for e in scenario.get('dist_events',[])})
    unsettled, final_errors, settling = 0, [], []
    for i, st in enumerate(starts):
        en = starts[i+1] if i+1 < len(starts) else t[-1]+dt
        mask = (t >= st) & (t < en)
        err = np.abs(sp[mask]-y[mask]); ts=t[mask]
        outside = np.flatnonzero(err > .05)
        val = 0. if not len(outside) else (np.inf if outside[-1]==len(err)-1 else ts[outside[-1]+1]-ts[0])
        unsettled += int(np.isinf(val)); settling.append(val)
        final_errors.append(float(err[-1]))
    return dict(IAE=float(np.sum(np.abs(sp-y))*dt), ISE=float(np.sum((sp-y)**2)*dt),
        max_abs_error=float(np.max(np.abs(sp-y))), actuator_TV=float(np.abs(np.diff(u)).sum()),
        saturation_fraction=float(np.mean((u<=1e-10)|(u>=FMAX-1e-10))),
        segments=len(starts), unsettled_segments=unsettled,
        worst_settling_s=float(max(settling)), worst_final_error=max(final_errors),
        pH_min=float(y.min()), pH_max=float(y.max()))

def simulate(job, keep_trace=False):
    scenario = (DEV if job['set']=='development' else TEST)[job['scenario']]
    dt = job.get('dt', 5.)
    n = max(800, int(scenario['tend']/20)) if dt == 'legacy' else int(round(scenario['tend']/dt))+1
    actual_dt = scenario['tend']/(n-1)
    delay = NS['theta']*job.get('delay_factor', 1.)
    SETTINGS.update(solver=job.get('solver','exact'), rtol=job.get('rtol',1e-9), atol=job.get('atol',1e-12))
    NS.update(NOMINAL)
    # Initial actuator bias and bounds are nominal, not recomputed from true chemistry.
    initial_flow = NOMINAL_FLOW(scenario['initial_pH'], 0.)
    for key in NOMINAL:
        NS[key] = NOMINAL[key]*job.get(key+'_factor',1.)
    NS['solve_Qb_for_pH'] = lambda *args, **kwargs: initial_flow
    NS['x_to_pH'] = stable_ph if job.get('stable',True) else LEGACY_PH
    rng = np.random.default_rng(job.get('seed',0))
    noise = rng.normal(0,job.get('noise_sd',0.),n)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res=NS['audited_sim'](**DESIGNS[job['design']], **scenario,
                n_points=n, theta_dead=delay, F_max=FMAX, noise=noise,
                delay_mode=job.get('delay_mode','legacy'), plot=False)
        row=dict(job, dt_actual=actual_dt, requested_delay=delay,
            implemented_delay=delay if job.get('delay_mode')=='fractional' else max(1,int(delay/actual_dt))*actual_dt,
            status='ok', **metrics(res,scenario))
        if keep_trace: return row,res
        return row
    except Exception as e:
        return dict(job,status='failed',failure=repr(e),IAE=np.nan)
    finally:
        NS.update(NOMINAL); NS['solve_Qb_for_pH']=NOMINAL_FLOW

def jobs():
    numerical=[]
    for design in DESIGNS:
        for scenario in DEV:
            base=dict(set='development',design=design,scenario=scenario)
            for label,options in [
                ('legacy_RK45',dict(dt='legacy',solver='rk',rtol=1e-3,atol=1e-6,stable=False)),
                ('tight_RK45',dict(dt='legacy',solver='rk',rtol=1e-9,atol=1e-12,stable=False)),
                ('exact_legacy',dict(dt='legacy',stable=False)),
                ('stable_legacy',dict(dt='legacy')),
                *[(f'sample_{dt}',dict(dt=dt)) for dt in [20.,10.,5.,2.5]],
                ('fractional_5',dict(dt=5.,delay_mode='fractional'))]:
                numerical.append(dict(base,condition=label,**options))
    conditions=[('nominal',{}),('delay_minus30',dict(delay_factor=.7)),
        ('delay_plus50',dict(delay_factor=1.5)),('acid_plus15',dict(Ca_factor=1.15)),
        ('acid_minus15',dict(Ca_factor=.85)),('base_minus15',dict(Cb_factor=.85)),
        ('base_plus15',dict(Cb_factor=1.15)),('Kw_half',dict(Kw_factor=.5)),
        ('Kw_double',dict(Kw_factor=2.)),('volume_plus20',dict(V_factor=1.2)),
        ('volume_minus20',dict(V_factor=.8))]
    for sd in [.01,.05]:
        for seed in [104729,130363,155921,196613,262147]:
            conditions.append((f'noise_{sd}_seed{seed}',dict(noise_sd=sd,seed=seed)))
    for seed in [104729,130363,155921,196613,262147]:
        conditions.append((f'combined_seed{seed}',dict(noise_sd=.05,seed=seed,
            delay_factor=1.5,Ca_factor=1.15,Cb_factor=.85,Kw_factor=2.,V_factor=1.2)))
    unseen=[dict(set='unseen',design=d,scenario=s,condition=c,dt=5.,**kw)
            for d in DESIGNS for s in TEST for c,kw in conditions]
    return numerical,unseen

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    numerical,unseen=jobs()
    protocol=dict(design_manifest_sha256=hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
        simulator_sha256=hashlib.sha256(SOURCE.encode()).hexdigest(),
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        theta=NS['theta'],F_max=FMAX,plant=NOMINAL,test_scenarios=TEST,
        numerical_jobs=numerical,unseen_jobs=unseen,
        criteria=dict(IAE_relative_change_percent=1.,trace_max_difference_pH=.01,
                      settling_band_pH=.05),
        note='Written before evaluation. Unseen means unused in historical tuning; no random population claim.')
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2))
    (OUT/'instrumented_original_simulator.py').write_text(AUDITED_SOURCE)
    # Verify zero-noise instrumentation matches unmodified original at defaults.
    check=dict(set='development',scenario='regulatory',design='reference_PI/E_original5_retuned',
               dt='legacy',solver='rk',rtol=1e-3,atol=1e-6,stable=False)
    row,aud=simulate(check,True)
    with contextlib.redirect_stdout(io.StringIO()):
        NS['x_to_pH']=LEGACY_PH
        original=NS['run_closed_loop_sim'](**DESIGNS[check['design']],**DEV['regulatory'],plot=False)
    for key in ['pH','F']:
        np.testing.assert_allclose(aud[key],original[key],rtol=1e-8,atol=1e-8)
    (OUT/'validation.txt').write_text('Zero-noise instrumentation versus unmodified original: PASS (rtol=atol=1e-8)\n')
    for name,batch in [('numerical_results',numerical),('unseen_results',unseen)]:
        rows=[]; start=time.time()
        with ProcessPoolExecutor(max_workers=4) as pool:
            for i,row in enumerate(pool.map(simulate,batch,chunksize=3)):
                rows.append(row)
                if (i+1)%30==0:
                    print(name,i+1,'/',len(batch),'elapsed',round(time.time()-start),flush=True)
                    pd.DataFrame(rows).to_csv(OUT/(name+'.partial.csv'),index=False)
        pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
        print(name,'finished',flush=True)

if __name__=='__main__': main()

"""Follow-up numerical checks; no controller tuning or test-scenario changes."""
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from run_validation import OUT, DESIGNS, DEV, TEST, simulate, jobs

def main():
    extra=[]
    for d in DESIGNS:
        for s in DEV:
            for dt in [1.25,.625]:
                extra.append(dict(set='development',design=d,scenario=s,condition=f'sample_{dt}',dt=dt))
            for dt in [2.5,1.25]:
                extra.append(dict(set='development',design=d,scenario=s,condition=f'fractional_{dt}',dt=dt,delay_mode='fractional'))
    # All unseen comparisons repeated at half the control interval, same frozen designs.
    _,unseen=jobs()
    refined=[dict(j,dt=2.5) for j in unseen]
    (OUT/'refinement_protocol.json').write_text(json.dumps(dict(numerical=extra,unseen=refined),indent=2))
    for name,batch in [('numerical_refinement',extra),('unseen_refinement',refined)]:
        rows=[]
        with ProcessPoolExecutor(max_workers=4) as pool:
            for i,row in enumerate(pool.map(simulate,batch,chunksize=3)):
                rows.append(row)
                if (i+1)%60==0: print(name,i+1,'/',len(batch),flush=True)
        pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
    # State trajectory agreement is more stringent than a whole-run integral.
    checks=[]
    for d in ['reference_PI/A_original3','reference_PI/E_original5_retuned','IMC_PID/E_original5_retuned']:
        for s in DEV:
            base=dict(set='development',design=d,scenario=s,dt='legacy',stable=False)
            _,legacy=simulate(dict(base,solver='rk',rtol=1e-3,atol=1e-6),True)
            _,tight=simulate(dict(base,solver='rk',rtol=1e-9,atol=1e-12),True)
            _,exact=simulate(dict(base,solver='exact'),True)
            for label,trace in [('default_to_tight',legacy),('exact_to_tight',exact)]:
                checks.append(dict(design=d,scenario=s,comparison=label,
                    max_pH_difference=float(np.max(np.abs(trace['pH']-tight['pH'])))))
    pd.DataFrame(checks).to_csv(OUT/'trajectory_checks.csv',index=False)

if __name__=='__main__':main()

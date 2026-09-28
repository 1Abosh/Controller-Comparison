"""Targeted checks of the new uncertainty and delay instrumentation."""
import hashlib
import json
import numpy as np
from run_validation import NS, OUT, SOURCE, MANIFEST, SETTINGS, advance, simulate

protocol=json.loads((OUT/'protocol.json').read_text())
assert protocol['design_manifest_sha256']==hashlib.sha256(MANIFEST.read_bytes()).hexdigest()
assert protocol['simulator_sha256']==hashlib.sha256(SOURCE.encode()).hexdigest()
base=dict(set='unseen',scenario='reverse_off_anchor',design='reference_PI/B_analytic3',dt=5.,
          delay_factor=100./NS['theta'])
_,a=simulate(dict(base,delay_mode='legacy'),True)
_,b=simulate(dict(base,delay_mode='fractional'),True)
np.testing.assert_allclose(a['pH'],b['pH'],rtol=1e-8,atol=1e-8)
_,a=simulate(dict(base,noise_sd=.05,seed=104729),True)
_,b=simulate(dict(base,noise_sd=.05,seed=104729),True)
np.testing.assert_array_equal(a['pH'],b['pH'])
max_error=0.
for x in [-.001,-1e-7,0.,1e-7,.001]:
    for q in [0.,5.,10.]:
        def fun(t,state,q,ca):
            return [(NS['Qa']*ca-q*NS['Cb']-(NS['Qa']+q)*state[0])/NS['V']]
        SETTINGS.update(solver='exact')
        a=advance(fun,[x],q,.01,20).y[0,-1]
        SETTINGS.update(solver='rk',rtol=1e-12,atol=1e-15)
        b=advance(fun,[x],q,.01,20).y[0,-1]
        max_error=max(max_error,abs(a-b))
assert max_error<1e-12
with (OUT/'validation.txt').open('a') as f:
    f.write('Frozen manifest and original simulator hashes unchanged: PASS\n')
    f.write('Integer-valued delay equivalence, legacy versus fractional: PASS\n')
    f.write('Seeded measurement-noise reproducibility: PASS\n')
    f.write(f'Exact interval propagation versus tight RK45, 15 states/flows: PASS (max state error {max_error:.3g})\n')
print('Implementation checks passed')

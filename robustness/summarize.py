"""Create paired, reproducible summaries without selecting or retuning designs."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/robustness'

def main():
    n=pd.concat([pd.read_csv(OUT/'numerical_results.csv'),pd.read_csv(OUT/'numerical_refinement.csv')],ignore_index=True)
    u=pd.read_csv(OUT/'unseen_results.csv')
    v=pd.read_csv(OUT/'unseen_refinement.csv')
    assert n.status.eq('ok').all() and u.status.eq('ok').all() and v.status.eq('ok').all()
    pairs=[('legacy_RK45','tight_RK45'),('tight_RK45','exact_legacy'),('exact_legacy','stable_legacy'),
        ('sample_20.0','sample_10.0'),('sample_10.0','sample_5.0'),('sample_5.0','sample_2.5'),
        ('sample_2.5','sample_1.25'),('sample_1.25','sample_0.625'),
        ('sample_5.0','fractional_5'),('fractional_5','fractional_2.5'),('fractional_2.5','fractional_1.25')]
    p=n.pivot(index=['design','scenario'],columns='condition',values='IAE')
    num=[]
    for a,b in pairs:
        change=(p[a]/p[b]-1).abs()*100
        num.append(dict(comparison=a+' to '+b,max_IAE_difference_percent=change.max(),
            median_IAE_difference_percent=change.median(),above_one_percent=int(sum(change>1)),n=len(change)))
    pd.DataFrame(num).to_csv(OUT/'numerical_summary.csv',index=False)
    nominal=u[u.condition=='nominal'].set_index(['design','scenario'])
    u['IAE_change_percent']=[100*(r.IAE/nominal.loc[(r.design,r.scenario),'IAE']-1) for r in u.itertuples()]
    baseline=u[u.design=='fixed/C1'].set_index(['scenario','condition']).IAE
    u['IAE_ratio_to_C1']=[r.IAE/baseline.loc[(r.scenario,r.condition)] for r in u.itertuples()]
    u.to_csv(OUT/'unseen_paired_results.csv',index=False)
    contrasts=[]
    for family in ['reference_PI','IMC_PID']:
        for a,b in [('A_original3','B_analytic3'),('B_analytic3','C_add7'),('C_add7','D_add8'),
                    ('D_add8','D2_resample_Ti_Td'),('D2_resample_Ti_Td','E_original5_retuned'),
                    ('D_add8','P_move_interior_Kc'),('A_original3','E_original5_retuned')]:
            left=u[u.design==family+'/'+a].set_index(['scenario','condition'])
            right=u[u.design==family+'/'+b].set_index(['scenario','condition'])
            for idx in left.index:
                contrasts.append(dict(family=family,before=a,after=b,scenario=idx[0],condition=idx[1],
                    IAE_change_percent=100*(right.loc[idx,'IAE']/left.loc[idx,'IAE']-1),
                    unsettled_before=int(left.loc[idx,'unsettled_segments']),
                    unsettled_after=int(right.loc[idx,'unsettled_segments'])))
    pd.DataFrame(contrasts).to_csv(OUT/'unseen_contrasts.csv',index=False)
    uv=u.merge(v,on=['design','scenario','condition'],suffixes=('_5s','_2p5s'))
    uv['IAE_absolute_change_percent']=100*(uv.IAE_5s/uv.IAE_2p5s-1).abs()
    uv[['design','scenario','condition','IAE_5s','IAE_2p5s','IAE_absolute_change_percent',
        'unsettled_segments_5s','unsettled_segments_2p5s']].to_csv(OUT/'unseen_sampling_comparison.csv',index=False)
    overall=u.groupby('design').agg(runs=('IAE','size'),mean_IAE=('IAE','mean'),worst_IAE=('IAE','max'),
        unsettled_segments=('unsettled_segments','sum'),total_segments=('segments','sum'),
        maximum_IAE_increase_percent=('IAE_change_percent','max'),max_actuator_TV=('actuator_TV','max'))
    overall.to_csv(OUT/'robustness_summary.csv')
    selected=['fixed/C1','reference_PI/A_original3','reference_PI/B_analytic3',
        'reference_PI/E_original5_retuned','IMC_PID/A_original3','IMC_PID/E_original5_retuned']
    print(pd.DataFrame(num).to_string(index=False))
    print('\nNOMINAL 5s\n',u[(u.condition=='nominal')&u.design.isin(selected)][['design','scenario','IAE','unsettled_segments','IAE_ratio_to_C1']].to_string(index=False))
    print('\nNOMINAL SAMPLING\n',uv[(uv.condition=='nominal')&uv.design.isin(selected)][['design','scenario','IAE_2p5s','IAE_absolute_change_percent','unsettled_segments_2p5s']].to_string(index=False))
    print('\nOVERALL\n',overall.loc[selected].to_string())
    print('\nTRAJECTORY\n',pd.read_csv(OUT/'trajectory_checks.csv').to_string(index=False))
    print('\nSAMPLING NOISELESS\n',uv[~uv.condition.str.startswith(('noise','combined'))].IAE_absolute_change_percent.describe().to_string())

if __name__=='__main__':main()

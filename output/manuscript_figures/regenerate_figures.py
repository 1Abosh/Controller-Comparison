"""Run with Python + numpy, pandas, scipy, matplotlib, IPython.
Model source is extracted unchanged from Lab_Scale_CSTR_Simulation_2_corrected.ipynb.
Selected tuning values are the saved notebook selections; no tuning sweep is rerun.
"""
from pathlib import Path
import contextlib, io, json, runpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent
plt.show=lambda: plt.close('all')
with contextlib.redirect_stdout(io.StringIO()):
    m=runpy.run_path(str(OUT/'model_from_notebook.py'),init_globals={'schedule_filter_tau':60.0})
from matplotlib import font_manager
# Register locally installed Office fonts when available; otherwise use system Calibri.
for font_path in Path('/Applications/Microsoft Word.app/Contents/Resources/DFonts').glob('*.ttf'):
    if font_path.name.lower().startswith('calibri'):
        font_manager.fontManager.addfont(str(font_path))
font_manager.findfont('Calibri', fallback_to_default=False)
plt.rcParams.update({'font.family':'Calibri','mathtext.fontset':'custom','mathtext.rm':'Calibri','mathtext.it':'Calibri:italic','mathtext.bf':'Calibri:bold','font.size':10,'axes.labelsize':11,'axes.titlesize':11,'legend.fontsize':9,'axes.linewidth':0.8,'xtick.direction':'out','ytick.direction':'out','xtick.major.size':4,'ytick.major.size':4,'svg.fonttype':'none','savefig.facecolor':'white','path.simplify':False})
params={'C1':.1,'C2':.02,'C2b':.1,'C3':5.,'C4':80.,'C4b':8.}
controllers={k:m['CONTROLLER_BUILDERS'][k][1](v) for k,v in params.items()}
colors=dict(zip(params,['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']))
styles=dict(zip(params,['-','--','-.','-', '--','-.']))
def finish(fig,name):
    for ax in fig.axes:
        ax.tick_params(bottom=True,left=True)
        for spine in ax.spines.values(): spine.set_visible(True)
        ax.grid(False, which='both', axis='both')
        ax.set_axisbelow(True)
    for ext in ['svg','png']: fig.savefig(OUT/f'{name}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
def val(x,p): return np.array([x.calculate_gain(z) for z in p]) if hasattr(x,'calculate_gain') else np.full_like(p,x)
p=np.linspace(4,10,1201)
g=np.array([m['local_gain_analytic'](z) for z in p])
# Figure 2: analytic steady-state relationship, evaluated densely through equivalence.
q=np.array([m['solve_Qb_for_pH'](z,5.) for z in np.linspace(2.01,11.7,2000)])
x=(m['Qa']*m['Ca']-q*m['Cb'])/(m['Qa']+q)
fig,axs=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
axs[0].plot(q,m['x_to_pH'](x),color='#0072B2');axs[0].set(xlabel=r'Base flow $Q_b$ (mL/s)',ylabel='Steady-state pH',title='(a) Titration curve')
axs[1].plot(p,g,color='#0072B2');axs[1].set(xlabel='Operating pH',ylabel=r'$K_p$ [pH/(mL/s)]',title='(b) Analytic local process gain')
finish(fig,'Figure_2_process_gain')
fig,axs=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
for ax,ids,title in zip(axs,[['C1','C2','C2b'],['C3','C4','C4b']],['(a) Reference-tuned PI','(b) IMC-derived controllers']):
    for cid in ids: ax.plot(p,val(controllers[cid]['Kc'],p),color=colors[cid],ls=styles[cid],label=cid)
    ax.set(xlabel='Operating pH',ylabel=r'$K_c$ [(mL/s)/pH]',title=title);ax.legend(frameon=False)
finish(fig,'Figure_3_controller_gains')
fig,axs=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
for ax,key,label,title in zip(axs,['Kc','Ti'],[r'$K_c$ [(mL/s)/pH]',r'$T_i$ (s)'],['(a) Proportional gain','(b) Integral time']):
    for cid in ['C1','C2','C2b']: ax.plot(p,val(controllers[cid][key],p),color=colors[cid],ls=styles[cid],label=cid)
    ax.set(xlabel='Operating pH',ylabel=label,title=title);ax.legend(frameon=False)
finish(fig,'Figure_4_PI_parameters')
traces={}; rows=[]
for cid,par in controllers.items():
    print('Simulating',cid,flush=True)
    with contextlib.redirect_stdout(io.StringIO()): r=m['run_closed_loop_sim'](**par,**m['standardised_scenario'],plot=False)
    traces[cid]=r
    pd.DataFrame({k:r[k] for k in ['t','SP','pH','F']}).to_csv(OUT/f'trace_{cid}.csv',index=False)
    grid=np.linspace(5,9,200);loop=val(par['Kc'],grid)*np.array([m['local_gain_analytic'](z) for z in grid])
    rows.append({'controller':cid,'selected_parameter':params[cid],'mean_step_IAE_standardised':m['mean_step_iae'](r),'loop_gain_max_min':loop.max()/loop.min()})
def responses(ids,name):
    fig,axs=plt.subplots(len(ids),2,figsize=(10,3*len(ids)),layout='constrained',squeeze=False)
    for row,cid in enumerate(ids):
        r=traces[cid]
        for col,ax in enumerate(axs[row]):
            ax.plot(r['t']/1000,r['pH'],lw=1.0,color=colors[cid],label=cid+' response')
            ax.step(r['t']/1000,r['SP'],where='post',ls='--',lw=1.5,color='black',label='Setpoint',zorder=3)
            for event in m['standardised_scenario']['dist_events']: ax.axvline(event['time']/1000,color='#888888',ls=':',lw=.8)
            ax.set(xlim=(0,170) if col==0 else (60,85),ylim=(min(r['pH'].min(),6)-.2,max(r['pH'].max(),8)+.2),xlabel='Time (10³ s)',ylabel='pH',title=f'({chr(97+row*2+col)}) {cid} — Standardised scenario\n'+('Full closed-loop response' if col==0 else 'Single-step response'))
            ax.legend(loc='lower right',frameon=True,framealpha=.95)
    finish(fig,name)
responses(['C2','C2b'],'Figure_5_GS_PI_responses')
responses(['C3','C4b'],'Figure_6_IMC_responses')
responses(['C4'],'Supplement_C4_response')
fig,axs=plt.subplots(1,2,figsize=(10,4.1),layout='constrained')
grid=np.linspace(5,9,200);kp=np.array([m['local_gain_analytic'](z) for z in grid])
for ax,ids,title in zip(axs,[['C1','C2','C2b'],['C3','C4','C4b']],['(a) Reference-tuned PI','(b) IMC-derived controllers']):
    for cid in ids: ax.plot(grid,val(controllers[cid]['Kc'],grid)*kp,color=colors[cid],ls=styles[cid],label=cid)
    ax.set(xlabel='Operating pH',ylabel=r'$K_c K_p$ (dimensionless)',title=title);ax.legend(frameon=False)
finish(fig,'Figure_7_loop_gain')
pd.DataFrame(rows).to_csv(OUT/'verification_metrics.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))

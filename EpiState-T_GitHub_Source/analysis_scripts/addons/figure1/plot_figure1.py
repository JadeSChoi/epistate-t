"""Rebuild the opening figure from verified Gill array values and published-case counts."""
from pathlib import Path
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
g=pd.read_csv(root/'working_example.csv')
d=pd.read_csv(root/'donor_summaries.csv')
colors=['#006A73','#AD5927']; markers=['o','D']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.4,'axes.spines.top':False,'axes.spines.right':False})
fig,(a,b)=plt.subplots(1,2,figsize=(8.45,3.55),gridspec_kw={'width_ratios':[1.12,1]});fig.subplots_adjust(left=.09,right=.985,top=.83,bottom=.26,wspace=.34)
baseline=g.iloc[0]['treated'];g=g.set_index('phase').loc[['before','during','after']]; a.axhspan(baseline-.05,baseline+.05,color='#E9F0F1');a.axhline(baseline,ls=':',lw=1,color='#7C9197')
for ix,arm in enumerate(['treated','mock']):
 v=g[arm].values;a.plot([0,1,2],v,'-' if ix==0 else '--',lw=2,marker=markers[ix],ms=6,c=colors[ix],mfc=colors[ix] if ix==0 else 'white',label=arm.capitalize(),zorder=4-ix)
a.set_xticks([0,1,2],['Before','During','After']);a.set_xlim(-.2,2.2);vv=g[['treated','mock']].values
a.set_ylim(max(0,min(vv.min(),baseline-.05)-.025),min(1,max(vv.max(),baseline+.05)+.035))
a.set_ylabel('Bulk methylation (beta)');a.set_title('A   cg24947637 in donor O2',fontweight='bold',loc='left');a.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,fontsize=9)
for ix,arm in enumerate(['treated','mock']):
 vals=d.loc[d.arm==arm].set_index('donor').loc[['O1','O2','O3'],'pct_all'].to_numpy()
 b.bar(np.arange(3)+(ix-.5)*.32,vals,width=.30,color=colors[ix] if ix==0 else 'white',edgecolor=colors[ix],linewidth=1.3,hatch='//' if ix else '',label=arm.capitalize())
b.set_xticks([0,1,2],['O1','O2','O3']);b.set_ylim(0,10);b.set_ylabel('Returned sites (% of all CpGs)');b.set_title('B   Matched cultures',fontweight='bold',loc='left');b.legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2,fontsize=9)
for x in a,b:
 x.spines['left'].set_color('#667B80');x.spines['bottom'].set_color('#667B80');x.tick_params(length=3)
for ext in ('png','pdf'): fig.savefig(root/f'Figure1.{ext}',dpi=325,bbox_inches='tight',pad_inches=.07)

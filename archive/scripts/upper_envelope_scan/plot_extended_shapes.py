from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
R=Path.cwd()/'runs/20261030_upper_envelope_scan';data=json.loads((R/'extended_shapes_plot_data.json').read_text());fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True);cmap=plt.get_cmap('viridis');norm=Normalize(11,12)
for ax,era in zip(axes.flat,['2022','2022EE','2023','2023BPix']):
    for q in data['rows']:
        if q['era']==era:ax.plot(data['x'],q['density'],color=cmap(norm(q['mass'])),lw=1.1,ls='--' if q['mass']==11 else '-')
    ax.set_title(era);ax.set_xlabel('Massa ricostruita [GeV]');ax.set_ylabel('Densità normalizzata in 6–12 [GeV⁻¹]');ax.grid(alpha=.2);ax.set_xlim(9,12)
fig.suptitle('DCB parametriche: ipotesi 11,0–12,0 GeV, passo 0,1\nTratteggio: riferimento a 11; curve oltre 11 estrapolate, nuisance nominali',fontsize=13)
fig.subplots_adjust(top=.88,right=.88,hspace=.3,wspace=.2);cax=fig.add_axes([.91,.15,.018,.65]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='Massa ipotizzata [GeV]')
fig.savefig(R/'extended_signal_shapes.png',dpi=160);fig.savefig(R/'extended_signal_shapes.svg');plt.close(fig)

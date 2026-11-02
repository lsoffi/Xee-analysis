"""Display saved, converged B-only C0 fits without performing new fits."""
from pathlib import Path
import json
import os
import numpy as np

P = Path(__file__).resolve().parents[2]
R = P / 'runs/20261102_c0_background_yields'
SOURCE = P / 'runs/20261031_separate_sr_c0_limits'
os.environ.setdefault('MPLCONFIGDIR', str(P / 'work/mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run():
    R.mkdir(parents=True, exist_ok=True)
    for side, title, degree in [('low', 'SR bassa + CR', 6), ('high', 'SR alta + CR', 4)]:
        source = SOURCE / side / 'm9.5'
        channels = json.loads((source / 'bins.json').read_text())['b']
        fit = json.loads((source / 'fit_results.json').read_text())['fit']['b']
        assert fit['status'] == 0 and fit['covQual'] == 3 and fit['continuity_passed']
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
        for i, (ax, era) in enumerate(zip(axes.flat, ['2022', '2022EE', '2023', '2023BPix']), 1):
            chunks = sorted([channels[f'SR_ch{i}'], channels[f'CR_ch{i}']], key=lambda q: q['mass'][0])
            x = np.concatenate([q['mass'] for q in chunks])
            data = np.concatenate([q['data'] for q in chunks])
            total = np.concatenate([q['total'] for q in chunks])
            ax.axvspan(9.093333333333334, 10.506666666666668, color='#EEEEEE', zorder=0)
            ax.errorbar(x, data, yerr=np.sqrt(data), fmt='.', color='#222222', ms=2.5, elinewidth=.6, label='Dati', zorder=4)
            ax.plot(x, total, color='#007E87', lw=1.8, label='Fit totale, r = 0', zorder=3)
            for prefix, color, label, style in [
                ('local_continuum_', '#555555', f'Continuo (SR grado {degree}, CR grado 4)', '--'),
                ('shape_upsilon1s_', '#D55E00', r'$\Upsilon(1S)$', '-'),
                ('shape_upsilon2s_', '#8A4CA8', r'$\Upsilon(2S)$', '-')]:
                y = np.concatenate([next(v for k, v in q['components'].items() if k.startswith(prefix)) for q in chunks])
                ax.plot(x, np.where(y > 0, y, np.nan), color=color, ls=style, lw=1.4, label=label)
            boundary = 9.093333333333334 if side == 'low' else 10.506666666666668
            ax.axvline(boundary, color='#888888', ls=':', lw=1)
            ax.set_title(era, loc='left', fontweight='bold')
            ax.set_yscale('log')
            ax.set_ylim(.1, max(total.max(), data.max()) * 2.5)
            ax.set_xlim(x[0] - .006666667, x[-1] + .006666667)
            ax.set_ylabel('Eventi / 13.33 MeV')
            ax.set_xlabel('Massa [GeV]')
            ax.grid(axis='y', alpha=.15, which='major')
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=10, bbox_to_anchor=(.5, .025))
        fig.suptitle(f'Fit senza segnale: {title} | quattro ere | raccordo C0\n'
                     r'$\mu$ delle due $\Upsilon$ e normalizzazioni libere; larghezze e code fissate', fontsize=14)
        fig.text(.5, .005, 'Area grigia: CR. Fit salvato riutilizzato: status = 0, covQual = 3. Nessun nuovo fit.', ha='center', fontsize=9)
        fig.tight_layout(rect=(0, .11, 1, .91))
        for ext in ['png', 'pdf']:
            fig.savefig(R / f'fit_bonly_{side}.{ext}', dpi=160)
        plt.close(fig)
    print(R)


if __name__ == '__main__':
    run()

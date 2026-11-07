"""Read-only comparison of existing degree 1/2 results with the degree 3 scan.

All generated artifacts live in degree3_comparison; this module never fits,
changes workspaces, or overwrites the historical degree 1/2 report.
"""
from pathlib import Path
import csv
import datetime
import html
import json
import math
import sys
import time

import numpy as np
import uproot
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P = Path.cwd()
R = P / 'runs/20261030_upper_envelope_scan'
O = R / 'degree3_comparison'
DEGREES = (1, 2, 3)
MASSES = tuple(k / 10 for k in range(105, 121))
ERAS = ('2022', '2022EE', '2023', '2023BPix')
COLORS = {1: '#D55E00', 2: '#0072B2', 3: '#009E73'}


def read_json(path, errors):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        errors.append({'path': str(path.relative_to(R)), 'reason': str(exc)})
        return None


def collect():
    errors, observed, bias, gof = [], [], [], []
    for p in sorted(R.glob('degree[123]/m*/results.json')):
        q = read_json(p, errors)
        if q is not None:
            observed.append(q)
    for p in sorted(R.glob('bias/degree[123]/m*/*/m*/result.json')):
        q = read_json(p, errors)
        if q is not None:
            bias.append(q)
    raw_gof_counts = {}
    for degree in DEGREES:
        paths = list((R / f'gof/degree{degree}').glob('toy*/m*/result.json'))
        raw_gof_counts[str(degree)] = len(paths)
    for q in observed:
        toys = []
        for p in sorted((R / f'gof/degree{q["degree"]}').glob(f'toy*/m{q["mass"]:.1f}/result.json')):
            t = read_json(p, errors)
            if t is not None:
                toys.append(t)
        if not toys:
            continue
        finite = [t for t in toys if math.isfinite(t['score'])]
        # Retain every finite recorded score, including failed-CR toys.
        n = len(finite)
        if not n:
            continue
        k = sum(t['score'] > q['gof_observed'] for t in finite)
        failed = sum(not t['CR_valid'] for t in finite)
        k_valid = sum(t['CR_valid'] and t['score'] > q['gof_observed'] for t in finite)
        lo = 0. if k == 0 else float(beta.ppf(.025, k, n-k+1))
        hi = 1. if k == n else float(beta.ppf(.975, k+1, n-k))
        gof.append(dict(degree=q['degree'], mass=q['mass'], observed=q['gof_observed'],
                        n=n, k=k, p=k/n, low=lo, high=hi, CR_invalid=failed,
                        p_failure_low=k_valid/n, p_failure_high=(k_valid+failed)/n,
                        r_boundaries=sum(t['r_boundary'] for t in finite),
                        unique_toy_ids=len({t['toy'] for t in finite}),
                        nonfinite=len(toys)-n))
    return observed, bias, gof, errors, raw_gof_counts


def save(fig, name):
    fig.tight_layout()
    fig.savefig(O / f'{name}.png', dpi=160)
    fig.savefig(O / f'{name}.svg')
    plt.close(fig)


def fit_source(q):
    return R / f'degree{q["degree"]}/m{q["mass"]:.1f}/joint/fitdiagnostics/fitDiagnostics_upper.root'


def fit_link(q, region):
    if q['degree'] == 3:
        return f'fits/degree3/m{q["mass"]:.1f}/fit_{region}.png'
    return f'../degree{q["degree"]}/m{q["mass"]:.1f}/fit_{region}.png'


def plot_degree3_fits(q):
    src = fit_source(q)
    if not src.exists():
        return
    D = O / f'fits/degree3/m{q["mass"]:.1f}'
    D.mkdir(parents=True, exist_ok=True)
    result_time = (R / f'degree3/m{q["mass"]:.1f}/results.json').stat().st_mtime
    if all((D/f'fit_{reg}.png').exists() and (D/f'fit_{reg}.png').stat().st_mtime > max(src.stat().st_mtime, result_time) for reg in ['SR', 'CR']):
        return
    with uproot.open(src) as f:
        for region in ['SR', 'CR']:
            fig = plt.figure(figsize=(13, 9))
            outer = fig.add_gridspec(2, 2, hspace=.33, wspace=.2)
            for i, (era, cell) in enumerate(zip(ERAS, outer), 1):
                d = f[f'shapes_fit_s/{region}_ch{i}']
                pred, edges = d['total'].to_numpy()
                bw = np.diff(edges)
                x = (edges[:-1]+edges[1:])/2
                data = d['data'].values('y') * bw
                pred = pred * bw
                gs = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
                ax, rx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
                ax.errorbar(x, data, yerr=np.sqrt(data), fmt='.', ms=2, color='black', lw=.5, label='Dati')
                ax.stairs(pred, edges, baseline=None, color=COLORS[3], lw=1.5, label='Totale S+B')
                ax.stairs(f[f'shapes_fit_b/{region}_ch{i}/total'].values()*bw, edges, baseline=None, color='gray', ls='--', label='Solo fondo')
                for proc, color, label in [('dy', '#0072B2', 'Envelope grado 3' if region == 'SR' else 'Chebyshev4'),
                                           ('upsilon1s', '#E69F00', 'Υ(1S)'), ('upsilon2s', '#CC79A7', 'Υ(2S)'), ('Zd', '#D55E00', 'Segnale DCB')]:
                    ax.stairs(d[proc].values()*bw, edges, baseline=None, color=color, label=label)
                ax.set_title(era)
                ax.set_ylabel('Eventi / 13,33 MeV')
                ax.tick_params(labelbottom=False)
                ax.legend(fontsize=7, ncol=2)
                ax.grid(alpha=.15)
                rx.axhline(0, color='gray', lw=.7)
                rx.plot(x, (data-pred)/np.sqrt(np.maximum(pred, 1e-10)), '.', color=COLORS[3], ms=2.5)
                rx.set_ylabel('(D−M)/√M')
                rx.set_xlabel('Massa [GeV]')
                rx.set_xlim(edges[0], edges[-1])
                rx.grid(alpha=.15)
            fig.suptitle(f'mX={q["mass"]:.1f} GeV | {region} | Envelope Cheb3/Bern3, CR4\nμ Υ dalla sola CR e fissate; DCB in SR e CR' + (' | Resa estrapolata' if q['mass'] > 11 else ''), fontsize=13)
            fig.subplots_adjust(top=.89, bottom=.07)
            fig.text(.5, .015, f'r̂={q["r"]:.4g}; Z locale={q["Z"]:.3f}; FD B/S status={q["fit"]["b"]["status"]}/{q["fit"]["s"]["status"]}, cov={q["fit"]["b"]["covQual"]}/{q["fit"]["s"]["covQual"]}' + (' | FLAG NUMERICO' if not q['numerical_ok'] else ''), ha='center', fontsize=9)
            for ext in ('png', 'svg'):
                fig.savefig(D/f'fit_{region}.{ext}', dpi=140)
            plt.close(fig)


def plot_fit_overlays(observed):
    lookup = {(q['degree'], q['mass']): q for q in observed}
    links = {}
    for mass in MASSES:
        if not all((degree, mass) in lookup for degree in DEGREES):
            continue
        qs = [lookup[degree, mass] for degree in DEGREES]
        sources = [fit_source(q) for q in qs]
        if not all(p.exists() for p in sources):
            continue
        D = O/'fit_overlays'; D.mkdir(exist_ok=True)
        target = D/f'm{mass:.1f}_SR.png'
        links[mass] = str(target.relative_to(O))
        if target.exists() and target.stat().st_mtime > max(p.stat().st_mtime for p in sources):
            continue
        files = [uproot.open(p) for p in sources]
        try:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            for i, (era, ax) in enumerate(zip(ERAS, axes.flat), 1):
                first = files[0][f'shapes_fit_s/SR_ch{i}']
                _, edges = first['total'].to_numpy(); bw = np.diff(edges); x = (edges[:-1]+edges[1:])/2
                data = first['data'].values('y')*bw
                ax.errorbar(x, data, yerr=np.sqrt(data), fmt='.', ms=2, color='black', lw=.5, label='Dati')
                for degree, f in zip(DEGREES, files):
                    d = f[f'shapes_fit_s/SR_ch{i}']
                    ax.stairs(d['total'].values()*bw, edges, baseline=None, color=COLORS[degree], lw=1.5, label=f'Grado {degree}: totale S+B')
                    ax.stairs(d['dy'].values()*bw, edges, baseline=None, color=COLORS[degree], ls='--', lw=1, label=f'Grado {degree}: continuo')
                ax.set_title(era);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 13,33 MeV');ax.grid(alpha=.15);ax.legend(fontsize=7)
            fig.suptitle(f'SR alta, mX={mass:.1f} GeV: confronto dei fit postfit\nLinee piene: totale S+B. Tratteggi: solo continuo. r e nuisance profilati per ogni grado.', fontsize=12)
            fig.tight_layout(rect=(0, 0, 1, .91));fig.savefig(target, dpi=150);fig.savefig(target.with_suffix('.svg'));plt.close(fig)
        finally:
            for f in files:f.close()
    return links


def plot_curves(observed, bias, gof):
    for name in ['limits', 'pvalues']:
        fig, ax = plt.subplots(figsize=(11, 5.8))
        for degree in DEGREES:
            rows = sorted((q for q in observed if q['degree'] == degree), key=lambda q:q['mass'])
            if not rows:continue
            x = [q['mass'] for q in rows]
            if name == 'limits':
                ax.plot(x, [q['limits']['observed'] for q in rows], 'o-', color=COLORS[degree], label=f'Grado {degree}: osservato')
                ax.plot(x, [q['limits']['expected_50'] for q in rows], '--', color=COLORS[degree], label=f'Grado {degree}: atteso')
                ax.fill_between(x, [q['limits']['expected_2p5'] for q in rows], [q['limits']['expected_97p5'] for q in rows], color=COLORS[degree], alpha=.055)
                ax.fill_between(x, [q['limits']['expected_16'] for q in rows], [q['limits']['expected_84'] for q in rows], color=COLORS[degree], alpha=.12)
            else:ax.plot(x, [q['p0'] for q in rows], 'o-', color=COLORS[degree], label=f'Envelope grado {degree}')
            bad = [q for q in rows if not q['numerical_ok']]
            if bad:ax.scatter([q['mass'] for q in bad], [q['limits']['observed'] if name == 'limits' else q['p0'] for q in bad], marker='x', color='black', s=80, zorder=5)
        ax.axvline(11, color='gray', ls=':');ax.axvspan(11, 12, color='#fff0d0', alpha=.35)
        ax.set_yscale('log');ax.set_xlim(10.48, 12.02);ax.set_xlabel('Massa ipotizzata [GeV]')
        ax.set_ylabel('Limite 95% CLs su r' if name == 'limits' else 'p-value locale (eccesso)')
        ax.set_title(f'SR alta + CR: confronto gradi 1, 2, 3 | {len(observed)}/48 configurazioni\n' + ('Bande native 1σ/2σ PROVVISORIE' if name == 'limits' else 'Significanze nominali; croci nere = flag numerico o minimo concorrente'))
        ax.grid(alpha=.2);ax.legend(fontsize=8, ncol=2 if name == 'limits' else 1);save(fig, name)
    if bias:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
        for ax, label, title in zip(axes, ['r0', 'rExpected', 'r10Expected'], ['r vero = 0', 'r vero = limite atteso', 'r vero = 10 × limite atteso']):
            for degree in DEGREES:
                rows = sorted((q for q in bias if q['degree'] == degree and q['label'] == label), key=lambda q:q['mass'])
                if not rows:continue
                ax.plot([q['mass'] for q in rows], [q['bias_sigma'] if q['valid'] else np.nan for q in rows], 'o-' if degree != 3 else 's--', color=COLORS[degree], label=f'Grado {degree}' + (' (chiusura)' if degree == 3 else ''))
                bad = [q for q in rows if not q['valid'] and q['bias_sigma'] is not None]
                ax.scatter([q['mass'] for q in bad], [q['bias_sigma'] for q in bad], marker='x', color=COLORS[degree], s=55)
            ax.axhspan(-.2, .2, color='gray', alpha=.15);ax.axhline(0, color='gray', lw=.6);ax.set_xlim(10.48, 12.02)
            ax.set_title(title);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('(r̂−r vero) / σHesse(r)');ax.grid(alpha=.2);ax.legend(fontsize=8)
        fig.suptitle(f'Asimov generati con Chebyshev3: {len(bias)}/144 fit\nGrado 3 = test di chiusura, NON bias da fondo alternativo. Croci = flag; μ ristimate in CR e poi fissate.', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .88));fig.savefig(O/'bias.png', dpi=160);fig.savefig(O/'bias.svg');plt.close(fig)
    if gof:
        fig, ax = plt.subplots(figsize=(11, 5.8))
        for degree in DEGREES:
            rows = sorted((q for q in gof if q['degree'] == degree), key=lambda q:q['mass'])
            if not rows:continue
            x=np.array([q['mass'] for q in rows]);p=np.array([q['p'] for q in rows]);lo=np.array([q['low'] for q in rows]);hi=np.array([q['high'] for q in rows])
            nmin=min(q['n'] for q in rows);nmax=max(q['n'] for q in rows);nlabel=str(nmin) if nmin == nmax else f'{nmin}–{nmax}'
            ax.errorbar(x+(degree-2)*.009, p, yerr=[p-lo, hi-p], fmt='o-', color=COLORS[degree], capsize=2, label=f'Grado {degree}: {nlabel} toy/punto')
        ax.axhline(.05, color='red', ls='--', label='0,05');ax.set_ylim(-.02, 1.02);ax.set_xlim(10.48, 12.02)
        ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('p GoF saturated, fit S+B con r libero')
        boundary_note=f'; r al bordo: fino a {max(q["r_boundaries"] for q in gof)} toy/punto' if any(q['r_boundaries'] for q in gof) else ''
        ax.set_title('GoF: generazione B-only del rispettivo modello; μ ristimate nella sola CR\nIntervalli MC 95%; toy correlati fra masse'+boundary_note,fontsize=11)
        ax.grid(alpha=.2);ax.legend(fontsize=9);save(fig, 'gof')


def export_csv(name, rows):
    if not rows:return
    with (O/name).open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def report():
    O.mkdir(parents=True, exist_ok=True)
    observed, bias, gof, errors, raw_gof_counts = collect()
    audit_path=R/'degree3/observed_audit.json'
    observed_audit=read_json(audit_path, errors) if audit_path.exists() else None
    audit_by_mass={q['mass']:q for q in observed_audit.get('results',[])} if observed_audit else {}
    observed_notes=[]
    for aq in audit_by_mass.values():
        if aq.get('refinement_sigma_unstable'):
            ref_path=R/f'degree3/m{aq["mass"]:.1f}/joint/refinement_degree3.json'
            ref=read_json(ref_path, errors) if ref_path.exists() else None
            trial_errors=[q['error'] for q in ref.get('trials',[]) if q.get('status')==0 and q.get('covQual')==3] if ref else []
            span=f'{min(trial_errors):.3f}–{max(trial_errors):.3f}' if trial_errors else 'variabile fra inizializzazioni'
            observed_notes.append(f'<li>A <b>{aq["mass"]:.1f} GeV</b>, σHesse(r)={span} fra soluzioni con likelihood coincidente: variazione relativa {100*aq["refinement_sigma_relative_spread"]:.1f}%. Il valore centrale e la significanza profilata sono riprodotti, ma <b>l’errore Hesse non è validato</b>. Questo avviso non viene trasformato automaticamente in un flag sul p-value.</li>')
        if aq.get('expected_warnings'):
            text=' '.join(aq['expected_warnings'])
            observed_notes.append(f'<li>A <b>{aq["mass"]:.1f} GeV</b>, il fit dell’Asimov usato per l’atteso non ritorna al valore r=0: <code>{html.escape(text)}</code>. Il limite atteso e le bande restano non validati; il crossing del limite osservato è verificato separatamente.</li>')
    for q in observed:
        if q['degree'] == 3:
            try:plot_degree3_fits(q)
            except (OSError, KeyError, ValueError) as exc:errors.append({'path':str(fit_source(q).relative_to(R)), 'reason':str(exc)})
    overlay_links = plot_fit_overlays(observed)
    plot_curves(observed, bias, gof)
    by_degree = {}
    for degree in DEGREES:
        a=[q for q in observed if q['degree']==degree];b=[q for q in bias if q['degree']==degree];g=[q for q in gof if q['degree']==degree]
        by_degree[str(degree)] = dict(observed=len(a), observed_target=16, observed_ok=sum(q['numerical_ok'] for q in a), bias=len(b), bias_target=48, bias_valid=sum(q['valid'] for q in b), gof_scores=sum(q['n'] for q in g), gof_target=1600, gof_complete_points=sum(q['n']==100 and q['unique_toy_ids']==100 for q in g))
    counts=dict(observed=len(observed), observed_target=48, bias=len(bias), bias_target=144, bias_valid=sum(q['valid'] for q in bias), gof_scores=sum(q['n'] for q in gof), gof_target=4800)
    methods = {}
    for method in ['fitdiagnostics', 'significance', 'limits', 'gof_observed']:
        receipts = [read_json(p, errors) for p in (R/'degree3').glob(f'm*/joint/{method}/combine.receipt.json')]
        methods[method] = dict(completed=sum(q is not None and q.get('exit_code') == 0 for q in receipts), failed=sum(q is not None and q.get('exit_code') != 0 for q in receipts), target=16)
    complete=all(q['observed']==16 and q['bias']==48 and q['gof_complete_points']==16 for q in by_degree.values()) and counts['gof_scores']==4800 and not errors
    now=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
    headline_numbers={}
    for degree in DEGREES:
        a=[q for q in observed if q['degree']==degree];g=[q for q in gof if q['degree']==degree]
        bs=[]
        for label in ['r0','rExpected','r10Expected']:
            samples=[q for q in bias if q['degree']==degree and q['label']==label];valid=[q for q in samples if q['valid']]
            bs.append(dict(injection=label,completed=len(samples),valid=len(valid),within_0p2_sigma=sum(abs(q['bias_sigma'])<.2 for q in valid),bias_sigma_min=min((q['bias_sigma'] for q in valid),default=None),bias_sigma_max=max((q['bias_sigma'] for q in valid),default=None)))
        zmax=max(a,key=lambda q:q['Z']) if a else None
        headline_numbers[str(degree)]=dict(maximum_local_Z=None if zmax is None else dict(mass=zmax['mass'],Z=zmax['Z'],p0=zmax['p0']),observed_limit_min=min((q['limits']['observed'] for q in a),default=None),observed_limit_max=max((q['limits']['observed'] for q in a),default=None),gof_p_min=min((q['p'] for q in g),default=None),gof_p_max=max((q['p'] for q in g),default=None),gof_failed_CR_per_mass_max=max((q['CR_invalid'] for q in g),default=0),gof_r_boundaries_total_scores=sum(q['r_boundaries'] for q in g),gof_points_with_r_boundary=[dict(mass=q['mass'],count=q['r_boundaries']) for q in g if q['r_boundaries']],bias=bs)
    summary=dict(updated_at=now, completed=complete, counts=counts, by_degree=by_degree, degree3_methods=methods, observed=observed, bias=bias, gof=gof, raw_gof_file_counts=raw_gof_counts, read_errors=errors,
                 headline_numbers=headline_numbers, degree3_observed_audit=observed_audit,
                 interpretation=dict(bias_generator='Chebyshev3 for all degrees', degree3_bias='Closure only: same degree/function family as generator; no independent alternative-background validation', injection_normalization='rExpected and r10Expected use the median expected limit of the SAME fitting degree; criterion matched, numerical rtrue differs between degrees', limit_bands='Native and provisional; upper quantile crossings not validated', extrapolation='Above 11 GeV original DCB formulas, logquadratic rates fitted 10–11; relative lnN at 11; no extrapolation nuisance', CR='Means estimated in CR-only B fit and frozen; normalization shared/floating; repeat first stage per generated sample', toy_reuse='100 B-only samples per degree, reused across 16 mass tests; 4800 scores are not 4800 independent generated samples'))
    (O/'summary.json').write_text(json.dumps(summary, indent=2))
    (O/'progress.json').write_text(json.dumps(dict(updated_at=now, completed=complete, counts=counts, by_degree=by_degree, degree3_methods=methods, read_errors=errors), indent=2))
    if complete:
        completion=dict(status='completed_with_documented_flags',completed_at=now,counts=counts,by_degree=by_degree,headline_numbers=headline_numbers,degree3_observed_audit_source='../degree3/observed_audit.json',degree3_sigma_hesse_unstable_masses=[q['mass'] for q in audit_by_mass.values() if q.get('refinement_sigma_unstable')],degree3_expected_warning_masses=[q['mass'] for q in audit_by_mass.values() if q.get('expected_warnings')],expected_bands='provisional',degree3_Asimov_interpretation='closure only, not alternative-background validation',injection_normalization='median expected limit of each fitting degree',input_results_modified=False)
        (O/'completion.json').write_text(json.dumps(completion,indent=2))
    export_csv('observed_results.csv', [dict(degree=q['degree'], mass=q['mass'], r=q['r'], r_error=q['r_error'], Z=q['Z'], p0=q['p0'], **q['limits'], numerical_ok=q['numerical_ok'], degree3_r_error_instability_documented=bool(audit_by_mass.get(q['mass'],{}).get('refinement_sigma_unstable')) if q['degree']==3 else False, degree3_expected_warning=' '.join(audit_by_mass.get(q['mass'],{}).get('expected_warnings',[])) if q['degree']==3 else '') for q in observed])
    export_csv('bias_results.csv', bias)
    export_csv('gof_results.csv', gof)
    progress_rows=''.join(f'<tr><td style="color:{COLORS[d]}">Grado {d}</td><td>{by_degree[str(d)]["observed"]}/16</td><td>{by_degree[str(d)]["bias"]}/48 ({by_degree[str(d)]["bias_valid"]} con qualità richiesta)</td><td>{by_degree[str(d)]["gof_scores"]}/1600</td></tr>' for d in DEGREES)
    fit_rows=[]
    for q in sorted(observed, key=lambda q:(q['mass'],q['degree'])):
        links=[]
        for reg in ['SR','CR']:
            link=fit_link(q,reg)
            links.append(f'<a href="{link}">{reg}</a>' if (O/link).exists() else f'{reg}: in preparazione')
        overlay=f'<a href="{overlay_links[q["mass"]]}">Sovrapposti 1/2/3</a>' if q['mass'] in overlay_links else 'In attesa del grado 3'
        quality='OK' if q['numerical_ok'] else 'FLAG'
        aq=audit_by_mass.get(q['mass'],{}) if q['degree']==3 else {}
        if aq.get('refinement_sigma_unstable'):quality+='; σHesse instabile'
        if aq.get('expected_warnings'):quality+=' osservati; atteso da verificare'
        fit_rows.append(f'<tr><td>{q["mass"]:.1f}</td><td>{q["degree"]}</td><td>{q["r"]:.4g}</td><td>{q["Z"]:.3f}</td><td>{q["p0"]:.3g}</td><td>{q["limits"]["observed"]:.4g}</td><td>{q["limits"]["expected_50"]:.4g}</td><td>{quality}</td><td>{" / ".join(links)}</td><td>{overlay}</td></tr>')
    bias_rows=[]
    for degree in DEGREES:
        for label in ['r0','rExpected','r10Expected']:
            rows=[q for q in bias if q['degree']==degree and q['label']==label];valid=[q for q in rows if q['valid']]
            bias_rows.append(f'<tr><td>{degree}{" (chiusura)" if degree==3 else ""}</td><td>{label}</td><td>{len(rows)}/16</td><td>{len(valid)}</td><td>{sum(abs(q["bias_sigma"])<.2 for q in valid)}</td></tr>')
    gof_rows=''.join(f'<tr><td>{q["degree"]}</td><td>{q["mass"]:.1f}</td><td>{q["k"]}/{q["n"]}</td><td>{q["p"]:.3f}</td><td>[{q["low"]:.3f}, {q["high"]:.3f}]</td><td>{q["CR_invalid"]}</td><td>{q["r_boundaries"]}</td><td>[{q["p_failure_low"]:.3f}, {q["p_failure_high"]:.3f}]</td></tr>' for q in gof)
    refresh='' if complete else '<meta http-equiv="refresh" content="60">'
    body=f'''<!doctype html><html lang="it"><meta charset="utf-8">{refresh}<title>SR alta: confronto gradi 1, 2, 3</title>
<style>body{{font:16px system-ui;max-width:1280px;margin:30px auto;padding:0 18px;color:#17324d}}img{{width:100%;height:auto}}table{{border-collapse:collapse;margin:20px 0}}th,td{{padding:8px;border:1px solid #ccd}}th{{background:#eef3f6}}.note{{padding:18px;background:#fff0d0;border-left:5px solid #b57900}}.status{{padding:16px;background:#eaf4f5}}a{{color:#00639a}}.scroll{{overflow-x:auto}}</style>
<h1>SR alta + CR: confronto dei gradi 1, 2 e 3</h1>
<p class="status">{'Calcoli completati con limitazioni documentate.' if complete else 'Confronto in aggiornamento: i risultati mancanti non sono interpolati.'}<br>Fit osservati {counts['observed']}/48; Asimov {counts['bias']}/144; GoF {counts['gof_scores']}/4800 score. Aggiornato {now}.</p>
<table><tr><th>Modello SR</th><th>Fit/limiti/p0</th><th>Asimov</th><th>GoF</th></tr>{progress_rows}</table>
<p>Nuovo grado 3, singoli metodi già completati: FitDiagnostics {methods['fitdiagnostics']['completed']}/16; significanza {methods['significance']['completed']}/16; limiti {methods['limits']['completed']}/16; GoF osservato {methods['gof_observed']['completed']}/16. Una configurazione entra nelle curve dopo l’estrazione di tutti i suoi metodi.</p>
<p>Stesse 16 masse, 10,5–12,0 GeV a passo 0,1; quattro ere combinate. SR 10,506667–12, CR 9,093333–10,506667 con Chebyshev4. Envelope Chebyshev/Bernstein dello stesso grado, con indice indipendente per era: 16 combinazioni. I risultati originali dei gradi 1/2 sono conservati e soltanto letti.</p>
<p class="note"><b>Grado 3: il test Asimov è una chiusura.</b> Il generatore rimane Chebyshev3, come nel confronto precedente. Fittare questi campioni con l’envelope Cheb3/Bern3 non costituisce un test di bias contro una famiglia di fondo alternativa. Un eventuale miglioramento non basta a dichiarare validato il grado 3. Nei due gradi inferiori il generatore contiene invece flessibilità aggiuntiva.</p>
<p>Le medie delle due Υ sono stimate nella sola CR B-only e poi congelate; le rese Υ sono condivise e libere nel simultaneo, σ/code MC fisse. Il segnale DCB è incluso in SR e CR, con unico r e nuisance originali. Il primo stadio CR è ripetuto in ciascun campione generato: il test riguarda l’intera procedura, compresa l’eventuale contaminazione della CR. Gli errori sui dati restano condizionati alle μ fissate.</p>
<h2>Limiti e p-value locali</h2><p class="note">Le bande attese native sono <b>provvisorie</b>: i quantili superiori non sono validati. I flag e i minimi concorrenti rimangono visibili. Le significanze sono locali, nominali e non corrette per lo scan; il GoF S+B non sostituisce il controllo di segnale spurio.</p>
{('<div class="note"><b>Controlli aggiuntivi del grado 3.</b><ul>'+''.join(observed_notes)+'</ul><p><a href="../degree3/observed_audit.json">Audit osservati del grado 3</a>. I controlli ordinari status/covarianza/positività/crossing CLs non equivalgono alla validazione degli errori o delle bande.</p></div>') if observed_notes else ''}
<img src="limits.png" alt="Limiti osservati e attesi per i tre gradi"><img src="pvalues.png" alt="P-value locali per i tre gradi">
<h2>Bias Asimov e chiusura</h2><p>Tre iniezioni: r=0, mediana del limite atteso del rispettivo modello, 10×mediana. <b>rExpected e r10Expected usano la mediana dello stesso grado che viene fittato</b>: stesso criterio di iniezione, ma r vero numericamente diverso fra i gradi 1/2/3; il valore effettivo è nel CSV. Non sono ensemble di bias e non misurano copertura. Croci: fit con flag; fascia ±0,2σ indicativa. Il denominatore usa l’errore Hesse condizionale alle μ fissate.</p>'''
    if (O/'bias.png').exists():body+='<img src="bias.png" alt="Bias Asimov; grado tre come chiusura">'
    body+=f'<table><tr><th>Grado</th><th>Iniezione</th><th>Completati</th><th>Con qualità richiesta</th><th>Entro ±0,2σ fra i validi</th></tr>{"".join(bias_rows)}</table><h2>Goodness of fit</h2><p>Generazione B-only del rispettivo modello e fit S+B con r libero; statistica Combine saturated. 100 campioni per grado, riutilizzati sulle 16 masse: correlazioni fra punti. p è la frazione strettamente sopra lo score osservato; intervalli binomiali Clopper–Pearson MC al 95%. I toy con CR non convergente restano nel denominatore; la colonna finale mostra l’intervallo se il loro score viene trattato come ignoto. La colonna r al bordo segnala separatamente i fit che raggiungono il range tecnico: sono conservati, senza allargamento automatico dei limiti. Non sono intervalli di copertura dell’analisi.</p>'
    if (O/'gof.png').exists():body+='<img src="gof.png" alt="GoF per i tre gradi">'
    body+=f'<div class="scroll"><table><tr><th>Grado</th><th>Massa</th><th>Superamenti/toy</th><th>p</th><th>MC 95%</th><th>Flag CR</th><th>r al bordo</th><th>Intervallo da flag CR</th></tr>{gof_rows}</table></div>'
    body+=f'<h2>Tutti i fit e le sovrapposizioni</h2><p>Ogni fit mostra quattro ere, componenti e residui. I confronti sovrappongono totale S+B e continuo dei tre gradi. FitDiagnostics ammette r negativo: il ramo della significanza fisica può essere diverso nei casi segnalati.</p><div class="scroll"><table><tr><th>Massa</th><th>Grado</th><th>r̂</th><th>Z</th><th>p0</th><th>Limite obs</th><th>Atteso</th><th>Numerica</th><th>Fit separati</th><th>Confronto SR</th></tr>{"".join(fit_rows)}</table></div>'
    body+='<h2>Segnale e tracciabilità</h2><p>Sopra 11 GeV: DCB dalle formule parametriche originali; rese estrapolate con log-quadratica 10–11 ancorata a 11. lnN relative mantenute a 11; nessuna nuisance aggiunta per l’estrapolazione. I valori nominali e le accettanze regionali sono gli stessi usati per i gradi 1/2. Fonti: directory degreeN, bias/degreeN e gof/degreeN del run originale.</p><p><a href="../signal_rate_extrapolation.png">Rese estrapolate</a> · <a href="../extended_signal_shapes.png">DCB sopra 11</a> · <a href="../index.html">Report storico gradi 1/2</a> · <a href="summary.json">Risultati JSON</a> · <a href="observed_results.csv">Fit/limiti/p0 CSV</a> · <a href="bias_results.csv">Asimov CSV</a> · <a href="gof_results.csv">GoF CSV</a></p>'
    for audit_candidate in [R/'degree3_validation_audit.json',R/'degree3/degree3_validation_audit.json']:
        if audit_candidate.exists():body+=f'<p><a href="../{audit_candidate.relative_to(R)}">Audit finale aggiuntivo del grado 3</a></p>';break
    if complete:body+='<p><a href="completion.json">Metadati di completamento e numeri riepilogativi</a></p>'
    if errors:body+='<p class="note">File temporaneamente non leggibili durante aggiornamento: '+html.escape(json.dumps(errors,ensure_ascii=False))+'</p>'
    body+='</html>'
    (O/'index.html').write_text(body)
    print(json.dumps(dict(completed=complete, counts=counts, degree3=by_degree['3'], degree3_methods=methods, read_errors=len(errors))), flush=True)
    return summary


if __name__ == '__main__':
    if '--watch' in sys.argv:
        while True:
            try:
                if report()['completed']:
                    break
            except Exception as exc:
                print(type(exc).__name__, str(exc), flush=True)
            time.sleep(30)
    else:
        report()

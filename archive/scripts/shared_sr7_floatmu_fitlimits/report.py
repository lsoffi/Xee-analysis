# -*- coding: utf-8 -*-
"""Read-only incremental fit/limit report: shared SR degree7, floating Y means."""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import html
import json
import math
import os
from zoneinfo import ZoneInfo

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'
os.environ.setdefault('MPLCONFIGDIR',str(P/'work/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MASSES=[i/10 for i in range(60,121)]
ERAS=['2022','2022EE','2023','2023BPix']
REGIONS=[('SRL','SR bassa'),('CR','CR'),('SRH','SR alta')]
MU_NAMES=[f'Ymean_{res}_{era}' for res in ['upsilon1s','upsilon2s'] for era in ERAS]
LIMIT_KEYS=['observed','expected_2p5','expected_16','expected_50','expected_84','expected_97p5']

def read(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return default

def finite(x):
    return isinstance(x, (float, int)) and math.isfinite(x)

def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False))
    temp.replace(path)

def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for q in rows for k in q))
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

def fmt(x, digits=4):
    return f'{x:.{digits}g}' if finite(x) else '—'

def fit_ok(q,stage=None):
    stages=[stage] if stage else ['b','s']
    return all(q.get('fit',{}).get(k,{}).get('numerical_ok',False)
               and finite(q['fit'][k].get('edm')) and 0<=q['fit'][k]['edm']<.01 for k in stages)


def point_ok(q,stage=None):
    stages=[stage] if stage else ['b','s']
    return all(q.get('fit',{}).get(k,{}).get('point_valid',False)
               and finite(q['fit'][k].get('edm')) and 0<=q['fit'][k]['edm']<.01 for k in stages)


def limit_ready(q):
    lim=q.get('limits') or {}
    return all(finite(lim.get(k)) for k in LIMIT_KEYS)


def fit_notes(q):
    notes=[]
    if q.get('fit') and not fit_ok(q):notes.append('minimo valido; covarianza FD da verificare' if point_ok(q) else 'FLAG numerico FD: vedere JSON')
    if q.get('limits') and not q.get('expected_bands_verified',q.get('expected_validated',False)):
        notes.append('bande attese provvisorie')
    if q.get('status')!='complete':notes.append('scan parziale')
    return '; '.join(notes)


def clean(value):
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):return [clean(v) for v in value]
    return value


def collect():
    rows=[]
    for p in sorted((R/'observed').glob('m*/results.json')):
        q=read(p)
        if q and q.get('mass') in MASSES:rows.append(q)
    return sorted(rows,key=lambda q:q['mass'])


def mu_table(rows):
    result=[]
    for q in rows:
        pre=q.get('mu_prefit') or {}
        params={}
        for stage in ['b','s']:
            raw=q.get('fit',{}).get(stage,{}).get('parameters',[])
            params[stage]={item['name']:item for item in raw} if isinstance(raw,list) else raw
        for name in MU_NAMES:
            a=params['b'].get(name,{});b=params['s'].get(name,{})
            p=pre.get(name)
            if isinstance(p,dict):p=p.get('value')
            if p is None and not a and not b:continue
            av,bv=a.get('value'),b.get('value')
            result.append(dict(mass=q['mass'],parameter=name,pre=p,
                post_B=av,error_B=a.get('error'),post_SB=bv,error_SB=b.get('error'),
                shift_SB_minus_pre_MeV=1000*(bv-p) if finite(bv) and finite(p) else None,
                shift_SB_minus_B_MeV=1000*(bv-av) if finite(bv) and finite(av) else None,
                low=b.get('low',a.get('low')),high=b.get('high',a.get('high')),
                B_numerical_ok=fit_ok(q,'b'),SB_numerical_ok=fit_ok(q,'s'),B_point_valid=point_ok(q,'b'),SB_point_valid=point_ok(q,'s'),
                B_boundary=name in (q.get('fit',{}).get('b',{}).get('mean_boundary',[]) or []),
                SB_boundary=name in (q.get('fit',{}).get('s',{}).get('mean_boundary',[]) or []),
                source=f'observed/m{q["mass"]:.1f}/results.json'))
    return result


def limit_plot(rows):
    if not any(q.get('limits') for q in rows):return
    lookup={q['mass']:q for q in rows}
    data=[dict(mass=m,limits=lookup.get(m,{}).get('limits'),
               fit_ok=fit_ok(lookup.get(m,{})),crossing=lookup.get(m,{}).get('observed_limit_check'),
               expected_validated=lookup.get(m,{}).get('expected_bands_verified',lookup.get(m,{}).get('expected_validated',False))) for m in MASSES]
    signature=hashlib.sha256(json.dumps(clean(data),sort_keys=True).encode()).hexdigest()
    if read(R/'limits_plot_cache.json',{}).get('signature')==signature and all((R/f'limits.{ext}').exists() for ext in ['png','pdf']):return
    def vals(key):
        return np.array([(lookup.get(m,{}).get('limits') or {}).get(key,np.nan) for m in MASSES],dtype=float)
    fig,ax=plt.subplots(figsize=(12,6.5))
    ax.fill_between(MASSES,vals('expected_2p5'),vals('expected_97p5'),color='#F6D743',label='Atteso ±2σ')
    ax.fill_between(MASSES,vals('expected_16'),vals('expected_84'),color='#69BE46',label='Atteso ±1σ')
    present=np.isfinite(vals('expected_50'))
    isolated=present & ~np.r_[False,present[:-1]] & ~np.r_[present[1:],False]
    for lo,hi,color,width in [('expected_2p5','expected_97p5','#F6D743',7),('expected_16','expected_84','#69BE46',4)]:
        good=isolated & np.isfinite(vals(lo)) & np.isfinite(vals(hi))
        ax.vlines(np.array(MASSES)[good],vals(lo)[good],vals(hi)[good],color=color,lw=width)
    ax.plot(MASSES,vals('expected_50'),'k--',lw=1.4,label='Atteso mediano')
    ax.plot(MASSES,vals('observed'),'o-',ms=3,color='#0072B2',label='Osservato')
    for q in rows:
        value=(q.get('limits') or {}).get('observed')
        if finite(value) and not (fit_ok(q) and q.get('observed_limit_check',{}).get('verified',False)):
            ax.plot(q['mass'],value,'x',color='#C00000',ms=8,mew=1.6)
    ax.plot([],[],'x',color='#C00000',ms=7,label='Flag numerici FD / crossing osservato')
    ax.axvspan(9.093333333333334,10.506666666666668,color='#888888',alpha=.1,label='Ipotesi nella CR')
    ax.axvspan(11,12,color='#D99A25',alpha=.09,label='Resa segnale estrapolata')
    ax.set_yscale('log');ax.set_xlim(5.97,12.03);ax.set_xlabel('Massa ipotizzata [GeV]')
    ax.set_ylabel('Limite superiore 95% CLs su r');ax.grid(alpha=.2)
    n=sum(limit_ready(q) for q in rows)
    fig.suptitle(f'Continuo comune SR: envelope Cheb7/Bern7; CR Cheb4 | {n}/61 limiti\nOtto medie Υ libere nel fit simultaneo; sigma e code fissate',fontsize=13)
    provisional=any(q.get('limits') and not q.get('expected_bands_verified',q.get('expected_validated',False)) for q in rows)
    note='Bande native provvisorie: estremi non verificati; valori conservati senza correzioni.' if provisional else 'Quantili nativi; nessuna interpolazione delle ipotesi mancanti.'
    fig.text(.5,.035,note,ha='center',fontsize=9,color='#9A2418' if provisional else 'black')
    ax.legend(fontsize=8,ncol=3,loc='best');fig.subplots_adjust(top=.84,bottom=.15)
    for ext in ['png','pdf']:
        target=R/f'limits.{ext}';temp=R/f'.limits_{os.getpid()}.{ext}';fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(R/'limits_plot_cache.json',dict(signature=signature))

def plot_full_spectrum(mass, q, curves, source, target):
    """Display the joint model; never bridge the SR polynomial through the CR."""
    image = target/'fit_full_log.png'
    if not all(f'{reg}_ch{i}' in curves['s'] for reg, _ in REGIONS for i in range(1, 5)):
        return None
    result_path = R/f'observed/m{mass:.1f}/results.json'
    result_time = result_path.stat().st_mtime if result_path.exists() else 0
    if read(target/'layout_version.json',{}).get('version')==2 and image.exists() and image.stat().st_mtime > max(source.stat().st_mtime, result_time):
        return str(image.relative_to(R))
    fig = plt.figure(figsize=(14, 10))
    outer = fig.add_gridspec(2, 2, hspace=.32, wspace=.2)
    negative_signal = False
    for i, (era, cell) in enumerate(zip(ERAS, outer), 1):
        sub = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
        ax, rx = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
        ymax = 1.
        for region, _ in REGIONS:
            channel = f'{region}_ch{i}'
            item = curves['s'][channel]
            x = np.asarray(item['mass'], dtype=float)
            data = np.asarray(item['data'], dtype=float)
            pred = np.asarray(item['total'], dtype=float)
            width = item['bin_width']; edges = np.r_[x-width/2, x[-1]+width/2]
            ymax = max(ymax, float(np.max(data+np.sqrt(np.maximum(data, 0)))))
            positive_data = data > 0
            ax.errorbar(x[positive_data], data[positive_data], yerr=np.sqrt(data[positive_data]), fmt='.', color='black', ms=1.8, lw=.4, label='Dati' if region=='SRL' else None)
            ax.stairs(np.where(pred>0, pred, np.nan), edges, baseline=None, color='#009E73', lw=1.6, label='Totale S+B' if region=='SRL' else None)
            bonly = curves.get('b', {}).get(channel)
            if bonly:
                btotal = np.asarray(bonly['total'], dtype=float)
                ax.stairs(np.where(btotal>0, btotal, np.nan), edges, baseline=None, color='gray', ls='--', lw=1, label='Solo fondo' if region=='SRL' else None)
            for name, raw in item['components'].items():
                values = np.asarray(raw, dtype=float)
                if 'upsilon1s' in name: color, title = '#E69F00', 'Υ(1S)'
                elif 'upsilon2s' in name: color, title = '#CC79A7', 'Υ(2S)'
                elif 'Zd' in name:
                    color, title = '#D55E00', 'Segnale DCB (parte positiva)'
                    negative_signal |= bool(np.any(values < 0))
                else:
                    color, title = ('#7B3294', 'Continuo CR4') if region=='CR' else ('#0072B2', 'Continuo comune SR')
                show_label = region=='SRL' or title=='Continuo CR4'
                ax.stairs(np.where(values>0, values, np.nan), edges, baseline=None, color=color, lw=1.1, label=title if show_label else None)
            residual = np.full_like(pred, np.nan); ok = np.isfinite(pred) & (pred > 0)
            residual[ok] = (data[ok]-pred[ok])/np.sqrt(pred[ok])
            rx.plot(x, residual, '.', color='#009E73', ms=2)
        for axis in [ax, rx]:
            axis.axvspan(9.093333333333334, 10.506666666666668, color='#999999', alpha=.13)
            axis.axvline(9.093333333333334, color='#999999', lw=.6)
            axis.axvline(10.506666666666668, color='#999999', lw=.6)
            axis.set_xlim(6, 12); axis.grid(alpha=.12)
        ax.set_title(era); ax.set_yscale('log'); ax.set_ylim(.15, ymax*1.9)
        ax.set_ylabel('Eventi / 13,33 MeV'); ax.tick_params(labelbottom=False)
        rx.axhline(0, color='gray', lw=.7); rx.set_ylabel('(D−M)/√M'); rx.set_xlabel('Massa [GeV]')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, fontsize=8.2, loc='upper center', bbox_to_anchor=(.5, .925), frameon=False)
    fig.suptitle(f'mX={mass:.1f} GeV | Fit simultaneo dell’intero spettro\nSR bassa e alta: stessa intensità Cheb7/Bern7; μ Υ libere; fascia grigia: CR Cheb4 indipendente', fontsize=13)
    fig.subplots_adjust(top=.85, bottom=.12)
    fit = q.get('fit', {})
    line = f'r̂={fmt(q.get("r"))}; FD B/S status={fit.get("b", {}).get("status", "—")}/{fit.get("s", {}).get("status", "—")}, cov={fit.get("b", {}).get("covQual", "—")}/{fit.get("s", {}).get("covQual", "—")}'
    if fit_notes(q): line += ' | '+fit_notes(q)
    fig.text(.5, .04, line, ha='center', fontsize=9)
    log_note = 'Scala log: i bin a zero restano nei residui; componenti non positive non visualizzate.'
    if negative_signal: log_note += ' La DCB stimata è negativa: vedere i plot regionali lineari.'
    fig.text(.5, .018, log_note, ha='center', fontsize=8)
    fig.savefig(image, dpi=145); fig.savefig(image.with_suffix('.svg')); plt.close(fig)
    return str(image.relative_to(R))

def plot_bins(mass, q, errors):
    source = R/f'observed/m{mass:.1f}/bins.json'
    curves = read(source)
    if not curves or not curves.get('s'):
        return None
    target = R/f'fits/m{mass:.1f}'
    target.mkdir(parents=True, exist_ok=True)
    links = []
    result_path = R/f'observed/m{mass:.1f}/results.json'
    result_time = result_path.stat().st_mtime if result_path.exists() else 0
    for region, label in REGIONS:
        channels = [f'{region}_ch{i}' for i in range(1, 5)]
        if not all(ch in curves['s'] for ch in channels):
            continue
        image = target/f'fit_{region}.png'
        links.append(str(image.relative_to(R)))
        if read(target/'layout_version.json',{}).get('version')==2 and image.exists() and image.stat().st_mtime > max(source.stat().st_mtime, result_time):
            continue
        fig = plt.figure(figsize=(13, 9))
        outer = fig.add_gridspec(2, 2, hspace=.32, wspace=.2)
        for era, channel, cell in zip(ERAS, channels, outer):
            item = curves['s'][channel]
            x = np.asarray(item['mass'], dtype=float)
            data = np.asarray(item['data'], dtype=float)
            pred = np.asarray(item['total'], dtype=float)
            width = item['bin_width']
            edges = np.r_[x-width/2, x[-1]+width/2]
            sub = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
            ax, rx = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
            ax.errorbar(x, data, yerr=np.sqrt(np.maximum(data, 0)), fmt='.', color='black', ms=2, lw=.5, label='Dati')
            ax.stairs(pred, edges, baseline=None, color='#009E73', lw=1.5, label='Totale S+B')
            bonly = curves.get('b', {}).get(channel)
            if bonly:
                ax.stairs(bonly['total'], edges, baseline=None, color='gray', ls='--', label='Solo fondo')
            for name, values in item['components'].items():
                if 'upsilon1s' in name: color, title = '#E69F00', 'Υ(1S)'
                elif 'upsilon2s' in name: color, title = '#CC79A7', 'Υ(2S)'
                elif 'Zd' in name: color, title = '#D55E00', 'Segnale DCB'
                else: color, title = '#0072B2', 'Continuo CR4' if region == 'CR' else 'Continuo comune SR'
                ax.stairs(values, edges, baseline=None, color=color, label=title)
            ax.set_title(era); ax.set_ylabel(f'Eventi / {width*1000:.2f} MeV'); ax.tick_params(labelbottom=False)
            ax.grid(alpha=.15); ax.legend(fontsize=7, ncol=2)
            residual = np.full_like(pred, np.nan)
            ok = np.isfinite(pred) & (pred > 0)
            residual[ok] = (data[ok]-pred[ok])/np.sqrt(pred[ok])
            rx.axhline(0, color='gray', lw=.7); rx.plot(x, residual, '.', color='#009E73', ms=2.5)
            rx.set_ylabel('(D−M)/√M'); rx.set_xlabel('Massa [GeV]'); rx.grid(alpha=.15)
            ax.set_xlim(edges[0], edges[-1]); rx.set_xlim(edges[0], edges[-1])
        fig.suptitle(f'mX={mass:.1f} GeV | {label} | Envelope comune Cheb7/Bern7, CR4\nOtto μ Υ libere nel fit simultaneo; sigma/code fisse; DCB in tutte e tre le regioni', fontsize=13)
        fig.subplots_adjust(top=.89, bottom=.11)
        fit = q.get('fit', {})
        annotation = f'r̂={fmt(q.get("r"))}'
        if fit:
            annotation += f'; FD B/S status={fit.get("b", {}).get("status", "—")}/{fit.get("s", {}).get("status", "—")}, cov={fit.get("b", {}).get("covQual", "—")}/{fit.get("s", {}).get("covQual", "—")}'
        if fit_notes(q): annotation += ' | '+fit_notes(q)
        fig.text(.5, .025, annotation, ha='center', fontsize=9)
        fig.savefig(image, dpi=140); fig.savefig(image.with_suffix('.svg')); plt.close(fig)
    full = plot_full_spectrum(mass, q, curves, source, target)
    if full:
        links.insert(0, full)
    write_json(target/'layout_version.json',dict(version=2))
    return links

def report():
    R.mkdir(parents=True,exist_ok=True)
    rows=collect();lookup={q['mass']:q for q in rows};means=mu_table(rows)
    limit_plot(rows)
    errors=[];gallery={}
    for q in rows:
        try:
            links=plot_bins(q['mass'],q,errors)
            if links:gallery[q['mass']]=links
        except (OSError,KeyError,ValueError) as exc:
            errors.append(dict(mass=q['mass'],reason=str(exc)))
    counts=dict(target=61,results=len(rows),fitdiagnostics=sum(all(k in q.get('fit',{}) for k in ['b','s']) for q in rows),
        fit_numerically_ok=sum(fit_ok(q) for q in rows),fit_B_ok=sum(fit_ok(q,'b') for q in rows),fit_SB_ok=sum(fit_ok(q,'s') for q in rows),
        limits=sum(limit_ready(q) for q in rows),observed_crossing_verified=sum(q.get('observed_limit_check',{}).get('verified',False) for q in rows),
        provisional_expected_bands=sum(limit_ready(q) and not q.get('expected_bands_verified',q.get('expected_validated',False)) for q in rows),
        complete_results=sum(q.get('status')=='complete' for q in rows),plot_points=len(gallery),
        complete_galleries=sum(len(v)==4 for v in gallery.values()),mu_prepost_entries=len(means),mu_pre_values=sum(finite(q['pre']) for q in means))
    counts.update(point_pairs_valid=sum(point_ok(q) for q in rows),point_B_valid=sum(point_ok(q,'b') for q in rows),point_SB_valid=sum(point_ok(q,'s') for q in rows))
    terminal={m:read(R/f'observed/m{m:.1f}/completion.json') for m in MASSES}
    failures={m:read(R/f'observed/m{m:.1f}/worker_failure.json') for m in MASSES}
    remaining=read(R/'remaining_runner.json',{})
    counts['terminal_mass_outcomes']=sum(terminal[m] is not None or failures[m] is not None for m in MASSES)
    counts['worker_failures']=sum(failures[m] is not None for m in MASSES)
    production_closed=remaining.get('status')=='completed' and counts['terminal_mass_outcomes']==61
    complete=counts['complete_results']==61 and counts['fitdiagnostics']==61 and counts['limits']==61 and counts['complete_galleries']==61 and production_closed
    watcher_done=production_closed and not errors
    now=datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds')
    summary=dict(updated_at=now,completed=complete,production_closed=production_closed,watcher_done=watcher_done,counts=counts,results=rows,mu_prepost=means,rendering_errors=errors,
        interpretation=dict(model='Common global Cheb7/Bern7 SR intensity and independent Cheb4 CR, with all eight Upsilon means floating in the joint fit.',
            mu_prefit='Initial values, not measurements with Gaussian prior errors; pre/post shifts are not pulls.',
            bands='Native expected quantiles retained; provisional when not verified.',
            numerical='Completion of recorded computations is separate from numerical/model validity. Source numerical flags are retained.'))
    write_json(R/'report_summary.json',clean(summary));write_json(R/'progress.json',dict(updated_at=now,completed=complete,production_closed=production_closed,watcher_done=watcher_done,remaining_runner_status=remaining.get('status'),renderer_parent_pid=os.getppid(),counts=counts,rendering_errors=errors))
    observed_csv=[]
    for q in rows:
        a=q.get('fit',{}).get('b',{});b=q.get('fit',{}).get('s',{})
        observed_csv.append(dict(mass=q['mass'],status=q.get('status'),r=q.get('r'),r_error=q.get('r_error'),
            B_status=a.get('status'),B_covQual=a.get('covQual'),B_EDM=a.get('edm'),SB_status=b.get('status'),SB_covQual=b.get('covQual'),SB_EDM=b.get('edm'),
            B_ok=fit_ok(q,'b'),SB_ok=fit_ok(q,'s'),B_point_valid=point_ok(q,'b'),SB_point_valid=point_ok(q,'s'),**(q.get('limits') or {}),flags=q.get('flags',[])))
    csv_write(R/'results.csv',observed_csv);csv_write(R/'mu_prepost.csv',means)
    fit_rows=[];mean_sections=[]
    for mass in MASSES:
        q=lookup.get(mass,{});lim=q.get('limits') or {};fit=q.get('fit',{})
        labels={'fit_full_log.png':'Intero spettro (log)','fit_SRL.png':'SR bassa','fit_CR.png':'CR','fit_SRH.png':'SR alta'}
        links=' · '.join(f'<a href="{p}">{labels[p.rsplit("/",1)[-1]]}</a>' for p in gallery.get(mass,[])) or 'In attesa'
        parts=[]
        for stage in ['b','s']:
            a=fit.get(stage,{})
            parts.append(f'{stage.upper()}: {fmt(a.get("status"))}/{fmt(a.get("covQual"))}; EDM {fmt(a.get("edm"),3)}')
        quality='OK numerico' if fit_ok(q) else 'Stima centrale valida; covarianza da verificare' if point_ok(q) else 'Flag numerico' if fit else 'In attesa'
        if failures[mass] is not None:quality='Worker fallito: vedere dettagli'
        detail=f'<a href="observed/m{mass:.1f}/results.json">JSON</a>' if q else '—'
        if failures[mass] is not None:detail+=f' · <a href="observed/m{mass:.1f}/worker_failure.json">Errore worker</a>'
        fit_rows.append(f'<tr><td>{mass:.1f}</td><td>{links}</td><td>{fmt(q.get("r"))}</td><td>{fmt(lim.get("observed"))}</td><td>{fmt(lim.get("expected_50"))}</td><td>{"<br>".join(parts)}</td><td>{quality}</td><td><a href="#mu{mass:.1f}">μ Υ</a> · {detail}</td></tr>')
        cell=[a for a in means if a['mass']==mass]
        if cell:
            mr=[]
            for a in cell:
                label=a['parameter'].replace('Ymean_upsilon1s_','Υ(1S) ').replace('Ymean_upsilon2s_','Υ(2S) ')
                boundary='B ' if a['B_boundary'] else '';boundary+='S+B' if a['SB_boundary'] else ''
                mr.append(f'<tr><td>{label}</td><td>{fmt(a["pre"],6)}</td><td>{fmt(a["post_B"],6)} ± {fmt(a["error_B"],3)}</td><td>{fmt(a["post_SB"],6)} ± {fmt(a["error_SB"],3)}</td><td>{fmt(a["shift_SB_minus_pre_MeV"],4)}</td><td>{fmt(a["shift_SB_minus_B_MeV"],4)}</td><td>{boundary or "—"}</td></tr>')
            body='<table><thead><tr><th>Media</th><th>Pre [GeV]</th><th>Post B [GeV]</th><th>Post S+B [GeV]</th><th>Δ(S+B−pre) [MeV]</th><th>Δ(S+B−B) [MeV]</th><th>Bordo segnalato</th></tr></thead><tbody>'+''.join(mr)+'</tbody></table>'
            mean_sections.append(f'<details id="mu{mass:.1f}"><summary>Massa {mass:.1f} GeV · {len(cell)} medie · {quality}</summary><div class="scroll">{body}</div></details>')
        else:mean_sections.append(f'<div id="mu{mass:.1f}" class="pending">Massa {mass:.1f}: valori pre/post in attesa.</div>')
    refresh='' if watcher_done else '<meta http-equiv="refresh" content="60">'
    limits_panel='<p class="pending">In attesa dei primi limiti nativi.</p>'
    if (R/'limits.png').exists():limits_panel='<a href="limits.png"><img src="limits.png" alt="Limiti con bande attese"></a><p><a href="limits.pdf">PDF dei limiti</a></p>'
    configlink='<a href="configuration.json">Configurazione</a> · ' if (R/'configuration.json').exists() else ''
    csvlinks=(' · <a href="results.csv">Risultati CSV</a>' if rows else '')+(' · <a href="mu_prepost.csv">Medie Υ CSV</a>' if means else '')
    page=f'''<!doctype html>
<html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">{refresh}
<title>Xee — grado 7, medie Υ libere: fit e limiti</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;color:#192639;max-width:1300px;margin:auto;padding:25px;background:#f4f6fa}}h1{{font-size:29px}}h2{{font-size:22px}}section,.card{{background:white;padding:20px;border:1px solid #d8e1ec;border-radius:12px;margin:16px 0}}.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:15px}}.count{{font-size:28px;font-weight:700}}.pending{{color:#68798d}}.notice{{border-left:4px solid #c58724;padding:12px;background:#fff6e7}}img{{max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:8px;text-align:left;border-bottom:1px solid #dde4ec}}th{{background:#edf2f8}}.scroll{{overflow:auto}}details{{margin:10px 0}}summary{{cursor:pointer;font-weight:600}}a{{color:#0758a2}}@media(max-width:750px){{.cards{{grid-template-columns:1fr}}body{{padding:12px}}}}</style>
</head><body><h1>Fit e limiti: SR comuni di grado 7, medie Υ libere</h1>
<p>{'Esecuzioni completate; i flag numerici rimangono visibili.' if complete else 'Esecuzioni concluse con risultati mancanti o errori: vedere le tabelle.' if production_closed else 'Campagna in corso: risultati aggiunti appena disponibili.'}<br>Aggiornato: {now}</p>
<div class="cards"><div class="card"><div class="count">{counts['fitdiagnostics']} / 61</div>FitDiagnostics B e S+B<br>{counts['point_pairs_valid']} coppie con stime centrali valide<br>{counts['fit_numerically_ok']} con anche covarianze valide<br>B: {counts['fit_B_ok']} · S+B: {counts['fit_SB_ok']}</div>
<div class="card"><div class="count">{counts['limits']} / 61</div>Limiti con sei quantili nativi<br>{counts['observed_crossing_verified']} crossing osservati verificati<br>{counts['provisional_expected_bands']} bande attese provvisorie</div>
<div class="card"><div class="count">{counts['complete_galleries']} / 61</div>Gallerie complete delle quattro ere<br>Spettro completo + tre regioni<br>{counts['mu_prepost_entries']} valori μ pre/post disponibili</div></div>
<section><h2>Configurazione di questo scan</h2>
<table><thead><tr><th>Regione</th><th>Range [GeV]</th><th>Continuo</th></tr></thead><tbody>
<tr><td>SR bassa</td><td>6–9,093333</td><td rowspan="2">Envelope Chebyshev7/Bernstein7: un’intensità globale comune alle due SR, coefficienti e famiglia condivisi per era; rese regionali legate dagli integrali</td></tr>
<tr><td>SR alta</td><td>10,506667–12</td></tr><tr><td>CR</td><td>9,093333–10,506667</td><td>Chebyshev4 indipendente</td></tr></tbody></table>
<p>61 ipotesi di massa da 6 a 12 GeV a passi di 0,1 GeV; quattro ere e binning originale. <b>Le otto medie di Υ(1S) e Υ(2S) sono libere nel fit simultaneo</b>, con normalizzazioni condivise fra regioni. Le sigma e le code restano fissate. Il segnale DCB compare in tutte e tre le regioni con lo stesso r. La penalità dell’envelope è contata una volta per era.</p>
<p class="notice">Sopra 11 GeV sono usate le shape DCB parametriche e le normalizzazioni estrapolate già definite. La convergenza numerica e il completamento dei file non costituiscono da soli validazione del modello.</p></section>
<section><h2>Limiti 95% CLs</h2><p>Bande native attese a 1 e 2 sigma. I quantili non verificati restano esplicitamente provvisori; i valori mancanti non vengono interpolati. Le croci segnalano i problemi numerici già registrati nei fit o nel crossing osservato.</p>{limits_panel}</section>
<section><h2>Tutti i fit e le componenti</h2><p>Ogni plot contiene le quattro ere, le componenti del fit e il confronto B/S+B. I residui sono diagnostici. Nella vista logaritmica le componenti non positive non sono disegnate; restano visibili nei plot regionali lineari.</p>
<div class="scroll"><table><thead><tr><th>Massa</th><th>Gallerie</th><th>r̂</th><th>Limite osservato</th><th>Atteso mediano</th><th>Status/covQual; EDM</th><th>Qualità</th><th>Dettagli</th></tr></thead><tbody>{''.join(fit_rows)}</tbody></table></div></section>
<section><h2>Medie Υ prima e dopo il fit</h2><p>“Pre” indica il valore iniziale, non una misura con prior gaussiano: gli spostamenti non sono pull in unità di sigma. Gli errori post-fit sono quelli restituiti dal fit e vanno letti insieme a covQual e ai flag numerici. Le tabelle mantengono distinti B e S+B.</p>{''.join(mean_sections)}</section>
<section><h2>Dati e avanzamento</h2><p>{configlink}<a href="progress.json">Avanzamento JSON</a> · <a href="report_summary.json">Riepilogo JSON</a>{csvlinks}</p><p>Esiti di massa conclusi: {counts['terminal_mass_outcomes']}/61 · worker falliti: {counts['worker_failures']}. Il monitor termina quando il runner delle 60 masse rimanenti è concluso e tutti i 61 esiti, incluso il pilota, sono documentati.</p><p>Errori di rendering: {len(errors)}. Il report legge i prodotti della nuova campagna senza modificarli.</p></section>
</body></html>'''
    tmp=R/'index.html.tmp';tmp.write_text(page,encoding='utf-8');tmp.replace(R/'index.html')
    print(json.dumps(dict(updated_at=now,completed=complete,counts=counts,rendering_errors=len(errors))),flush=True)
    return complete


if __name__=='__main__':
    report()

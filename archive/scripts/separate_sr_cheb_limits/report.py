# -*- coding: utf-8 -*-
# Read-only native-limit report for two separate SR+CR analyses.
from pathlib import Path
import argparse,csv,datetime,hashlib,html,json,math,os
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_separate_sr_cheb_limits'
os.environ.setdefault('MPLCONFIGDIR',str(P/'work/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
KEYS=['observed','expected_2p5','expected_16','expected_50','expected_84','expected_97p5']
REGIONS=['low','high']
DEGREES={'low':6,'high':4}
LABELS={'low':'SR bassa','high':'SR alta'}
OUTER=[('low',n/10) for n in range(60,91)]+[('high',n/10) for n in range(106,121)]


def read(path,default=None):
    try:return json.loads(path.read_text())
    except (OSError,json.JSONDecodeError):return default


def finite(value):return isinstance(value,(int,float)) and math.isfinite(value)
def fmt(value):return f'{value:.5g}' if finite(value) else '—'


def clean(value):
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):return [clean(v) for v in value]
    return value


def write_json(path,value):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(clean(value),indent=2),encoding='utf-8');temp.replace(path)


def csv_write(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def load_manifest(path):
    manifest=read(path)
    if not isinstance(manifest,dict):return None,OUTER
    points=manifest.get('points',[])
    selected=[]
    for item in points:
        region=item.get('region');mass=item.get('mass')
        if region not in REGIONS or not finite(mass):raise ValueError('Invalid manifest point')
        pair=(region,round(float(mass),1))
        if abs(pair[1]-mass)>1e-8:raise ValueError('Manifest mass is not at a 0.1 GeV step')
        if pair in selected:raise ValueError('Duplicate manifest point')
        selected.append(pair)
    if not selected:raise ValueError('Manifest points list is empty')
    return manifest,selected


def collect():
    rows=[]
    for region in REGIONS:
        for path in sorted((R/region).glob('m*/results.json')):
            q=read(path)
            if not q:continue
            mass=q.get('mass',float(path.parent.name[1:]))
            if q.get('region',region)!=region:raise ValueError('Result region disagrees with its directory: '+str(path))
            lim=q.get('limits') or {};receipt=q.get('receipt') or {};check=q.get('observed_limit_check') or {}
            rows.append(dict(region=region,mass=float(mass),degree=q.get('degree',DEGREES[region]),
                **{key:lim.get(key) for key in KEYS},status=q.get('status'),
                six_quantiles=all(finite(lim.get(key)) for key in KEYS),
                method_complete=q.get('status')=='complete' and receipt.get('exit_code')==0 and all(finite(lim.get(key)) for key in KEYS),
                observed_crossing_verified=check.get('verified') is True,
                observed_curve_included=check.get('verified') is True and receipt.get('exit_code')==0,
                observed_display_status='crossing_verified' if check.get('verified') is True and receipt.get('exit_code')==0 else 'unvalidated_native_value',
                observed_crossing_CLs=check.get('CLs'),
                expected_bands_verified=q.get('expected_bands_verified',q.get('expected_validated',False)) is True,
                exit_code=receipt.get('exit_code'),flags=q.get('flags',[]),
                source=str(path.relative_to(R))))
    return sorted(rows,key=lambda q:(REGIONS.index(q['region']),q['mass']))


def plot(rows,points):
    if not any(finite(q['observed']) or finite(q['expected_50']) for q in rows):return
    signature=hashlib.sha256(json.dumps(clean(dict(rows=rows,points=points,plot_version=4)),sort_keys=True).encode()).hexdigest()
    if read(R/'limits_plot_cache.json',{}).get('signature')==signature and all((R/f'limits.{ext}').exists() for ext in ['png','pdf']):return
    fig,axes=plt.subplots(1,2,figsize=(14,6.8))
    for region,ax,color in zip(REGIONS,axes,['#0072B2','#7B3294']):
        region_rows=[q for q in rows if q['region']==region]
        lookup={q['mass']:q for q in region_rows}
        requested=[m for reg,m in points if reg==region]
        if not requested:
            ax.text(.5,.5,'Nessuna ipotesi richiesta',ha='center',transform=ax.transAxes);continue
        grid=[n/10 for n in range(round(min(requested)*10),round(max(requested)*10)+1)]
        def vals(key):return np.array([lookup.get(m,{}).get(key,np.nan) if m in requested else np.nan for m in grid],dtype=float)
        ax.fill_between(grid,vals('expected_2p5'),vals('expected_97p5'),color='#F6D743',label='Atteso ±2σ')
        ax.fill_between(grid,vals('expected_16'),vals('expected_84'),color='#69BE46',label='Atteso ±1σ')
        present=np.isfinite(vals('expected_50'));isolated=present & ~np.r_[False,present[:-1]] & ~np.r_[present[1:],False]
        for lo,hi,col,width in [('expected_2p5','expected_97p5','#F6D743',7),('expected_16','expected_84','#69BE46',4)]:
            good=isolated & np.isfinite(vals(lo)) & np.isfinite(vals(hi))
            ax.vlines(np.array(grid)[good],vals(lo)[good],vals(hi)[good],color=col,lw=width)
        ax.plot(grid,vals('expected_50'),'k--',marker='_',ms=7,lw=1.4,label='Atteso mediano')
        observed_curve=np.array([lookup.get(m,{}).get('observed',np.nan) if m in requested and lookup.get(m,{}).get('observed_curve_included',False) else np.nan for m in grid],dtype=float)
        ax.plot(grid,observed_curve,'o-',color=color,ms=3,label='Osservato: crossing verificato')
        bad=[q for q in region_rows if q['mass'] in requested and finite(q['observed']) and not q['observed_curve_included']]
        if bad:ax.plot([q['mass'] for q in bad],[q['observed'] for q in bad],'x',color='#C00000',ms=8,mew=1.7,zorder=5,label='Valore nativo non validato (×)')
        if any(9.093333333333334<=m<=10.506666666666668 for m in requested):
            ax.axvspan(9.093333333333334,10.506666666666668,color='#888888',alpha=.1,label='Ipotesi nella CR')
        n=sum(q['six_quantiles'] and q['mass'] in requested for q in region_rows)
        ax.set_title(f'{LABELS[region]} Chebyshev{DEGREES[region]} + CR Chebyshev4 | {n}/{len(requested)}')
        ax.set_yscale('log');ax.set_xlim(min(grid)-.06,max(grid)+.06);ax.set_xlabel('Massa ipotizzata [GeV]')
        ax.set_ylabel('Limite superiore 95% CLs su r');ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best')
    fig.suptitle('Quattro ere combinate: limiti 95% CLs su r da due analisi SR+CR separate\nOtto medie Υ libere; sigma/code fisse; segnale DCB nelle regioni selezionate',fontsize=13)
    provisional=any(q['six_quantiles'] and not q['expected_bands_verified'] for q in rows)
    note='Bande native provvisorie/non validate: valori conservati senza correzioni.' if provisional else 'Quantili nativi conservati; nessuna interpolazione delle masse mancanti.'
    fig.text(.5,.035,note+'\n× rosse: crossing osservato non verificato / errore; curva interrotta in quei punti.\nLa stessa CR compare nelle due analisi separate: queste curve non sono un limite combinato.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.81,bottom=.22,wspace=.24)
    for ext in ['png','pdf']:
        target=R/f'limits.{ext}';temp=R/f'.limits_{os.getpid()}.{ext}';fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(R/'limits_plot_cache.json',dict(signature=signature))


def report(manifest_path):
    R.mkdir(parents=True,exist_ok=True)
    manifest,points=load_manifest(manifest_path);rows=collect();lookup={(q['region'],q['mass']):q for q in rows}
    selected=[q for q in rows if (q['region'],q['mass']) in points]
    extras=[q for q in rows if (q['region'],q['mass']) not in points]
    plot(selected,points)
    counts={region:dict(target=sum(reg==region for reg,m in points),results=sum(q['region']==region for q in selected),limits=sum(q['region']==region and q['six_quantiles'] for q in selected),verified_crossings=sum(q['region']==region and q['observed_crossing_verified'] for q in selected),completed_methods=sum(q['region']==region and q['method_complete'] for q in selected),unvalidated_observed=sum(q['region']==region and finite(q['observed']) and not q['observed_curve_included'] for q in selected)) for region in REGIONS}
    all_available=bool(manifest) and not manifest.get('CR_masses_pending',False) and all(lookup.get(point,{}).get('method_complete',False) for point in points)
    completed=all_available
    remaining=read(R/'remaining_runner.json',{})
    production_closed=remaining.get('status')=='completed'
    stopped=(R/'stop_record.json').exists()
    watcher_done=completed or production_closed or stopped
    now=datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds')
    out=dict(updated_at=now,completed=completed,production_closed=production_closed,stopped_by_user=stopped,watcher_done=watcher_done,renderer_parent_pid=os.getppid(),manifest_available=manifest is not None,manifest_path=str(manifest_path),target=len(points),counts=counts,points=[dict(region=r,mass=m) for r,m in points],results=selected,unrequested_results=extras,
        interpretation='Separate low-SR+CR and high-SR+CR analyses, statistically correlated through their shared CR data. The common CR is not multiplied or combined. Overlapping CR mass hypotheses are shown for both declared likelihoods, without selecting a curve based on its observed limit. Native limits and bands are unchanged; missing masses are not interpolated.')
    write_json(R/'limits.json',out);csv_write(R/'limits.csv',selected)
    write_json(R/'progress.json',{k:out[k] for k in ['updated_at','completed','production_closed','stopped_by_user','watcher_done','renderer_parent_pid','manifest_available','target','counts']})
    tables=[]
    for region in REGIONS:
        cells=[]
        for reg,m in points:
            if reg!=region:continue
            q=lookup.get((reg,m),{})
            flags=html.escape('; '.join(q.get('flags',[]))) or '—'
            source=f'<a href="{q["source"]}">JSON</a>' if q else 'Non disponibile' if production_closed or stopped else 'In attesa'
            crossing='Verificato' if q.get('observed_curve_included') else 'NON VERIFICATO: valore non validato' if q else '—'
            row_style=' class="unvalidated"' if q and not q.get('observed_curve_included') else ''
            method_label=html.escape(str(q.get('status') or '—'))+'; exit '+str(q.get('exit_code','—'))
            cells.append(f'<tr{row_style}><td>{m:.1f}</td><td>{fmt(q.get("observed"))}</td><td>{fmt(q.get("expected_2p5"))}</td><td>{fmt(q.get("expected_16"))}</td><td>{fmt(q.get("expected_50"))}</td><td>{fmt(q.get("expected_84"))}</td><td>{fmt(q.get("expected_97p5"))}</td><td>{crossing}</td><td>{method_label}</td><td>{flags}</td><td>{source}</td></tr>')
        tables.append(f'<section><h2>{LABELS[region]} + CR</h2><div class="scroll"><table><thead><tr><th>Massa</th><th>Osservato</th><th>Atteso −2σ</th><th>−1σ</th><th>Mediano</th><th>+1σ</th><th>+2σ</th><th>Crossing obs.</th><th>Esito metodo</th><th>Flag sorgente</th><th>Risultato</th></tr></thead><tbody>'+''.join(cells)+'</tbody></table></div></section>')
    refresh='' if watcher_done else '<meta http-equiv="refresh" content="60">'
    plot_html='<p class="pending">In attesa dei primi limiti nativi.</p>'
    if (R/'limits.png').exists():plot_html='<a href="limits.png"><img src="limits.png" alt="Limiti separati SR bassa e alta con bande"></a><p><a href="limits.pdf">PDF del grafico</a></p>'
    manifest_note='Elenco dei punti definito dal manifest.' if manifest else 'Manifest ancora in attesa: piano esterno mostrato provvisoriamente, 31 masse 6–9 GeV e 15 masse 10,6–12 GeV. Nessun completamento dichiarato prima del manifest.'
    links='<a href="limits.json">Dati e configurazione JSON</a> · <a href="progress.json">Avanzamento</a>'
    if selected:links+=' · <a href="limits.csv">Valori CSV</a>'
    if (R/'completion_audit.json').exists():links+=' · <a href="completion_audit.json">Audit di completezza e crossing</a>'
    if manifest_path.exists():links+=f' · <a href="{os.path.relpath(manifest_path,R)}">Manifest</a>'
    page=f'''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">{refresh}<title>Xee — limiti separati SR bassa/alta</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1300px;margin:auto;padding:25px;background:#f4f6fa;color:#192639}}section,.card{{background:white;border:1px solid #d8e1ec;border-radius:12px;padding:20px;margin:16px 0}}.cards{{display:grid;grid-template-columns:1fr 1fr;gap:15px}}.count{{font-size:28px;font-weight:700}}.pending{{color:#68798d}}.notice{{background:#fff6e7;padding:12px;border-left:4px solid #c58724}}img{{max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:8px;text-align:left;border-bottom:1px solid #dde4ec}}th{{background:#edf2f8}}.scroll{{overflow:auto}}.unvalidated{{background:#fff0ec}}a{{color:#0758a2}}@media(max-width:700px){{.cards{{grid-template-columns:1fr}}body{{padding:12px}}}}</style></head><body>
<h1>Limiti separati: SR bassa + CR e SR alta + CR</h1><p>{'Tutti gli esiti di calcolo previsti sono disponibili; i valori con crossing non verificato non sono limiti validati.' if completed else 'Produzione interrotta: sono mostrati soltanto i risultati disponibili.' if stopped else 'Esecuzioni concluse con risultati mancanti o errori: vedere le tabelle.' if production_closed else 'Risultati in corso di produzione.'}<br>Aggiornato: {now}</p>
<div class="cards"><div class="card"><div class="count">{counts['low']['limits']} / {counts['low']['target']}</div>SR bassa: Chebyshev6 + CR Chebyshev4<br>{counts['low']['verified_crossings']} crossing osservati verificati · {counts['low']['unvalidated_observed']} valori non validati</div><div class="card"><div class="count">{counts['high']['limits']} / {counts['high']['target']}</div>SR alta: Chebyshev4 + CR Chebyshev4<br>{counts['high']['verified_crossings']} crossing osservati verificati · {counts['high']['unvalidated_observed']} valori non validati</div></div>
<section><h2>Configurazione</h2><p>SR bassa: 6–9,093333 GeV. SR alta: 10,506667–12 GeV. CR comune come definizione: 9,093333–10,506667 GeV. Ogni curva deriva da una propria likelihood SR+CR: <b>la CR non viene moltiplicata e le due curve non sono combinate</b>. Le due analisi usano gli stessi dati della CR e sono quindi statisticamente correlate.</p><p>Singole funzioni Chebyshev, grado 6 nella SR bassa e grado 4 nella SR alta; grado 4 nella CR in entrambe le analisi. Otto medie Υ libere, sigma e code fissate; segnale DCB in tutte le regioni selezionate. {manifest_note} Le ipotesi 9,1–10,5 GeV nella CR sono valutate in entrambe le analisi, mantenendo separate le due curve e le rispettive likelihood.</p><p class="notice">Quattro ere combinate in ciascuna likelihood. Quantili e bande nativi conservati; gli eventuali flag di ricerca del bracket restano visibili nelle tabelle. Le bande senza verifica documentata restano provvisorie/non validate. Sopra 11 GeV sono usate le DCB parametriche e le normalizzazioni di segnale estrapolate nello studio già autorizzato, senza aggiungere una nuova incertezza di estrapolazione. Il completamento dei calcoli non equivale a validazione del modello.</p></section>
<section><h2>Limiti 95% CLs su r</h2><p>Osservato, atteso mediano e bande attese a 1 e 2 sigma. Le croci rosse riportano i valori nativi con crossing osservato non verificato o errore di esecuzione: non sono limiti validati. La curva osservata si interrompe in quei punti, mentre CSV e tabelle conservano integralmente i numeri e i flag originali. I punti non ancora pronti restano mancanti. I pannelli usano la stessa unità r, con scale verticali logaritmiche proprie.</p>{plot_html}</section>
{''.join(tables)}<section><h2>Dati</h2><p>{links}</p><p>Risultati presenti fuori dall’elenco corrente: {len(extras)}. Il report legge i file sorgente senza modificarli.</p></section></body></html>'''
    temp=R/'index.html.tmp';temp.write_text(page,encoding='utf-8');temp.replace(R/'index.html')
    print(json.dumps(dict(updated_at=now,completed=completed,manifest_available=manifest is not None,target=len(points),counts=counts)),flush=True)
    return completed


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,default=R/'manifest.json');args=parser.parse_args();report(args.manifest)

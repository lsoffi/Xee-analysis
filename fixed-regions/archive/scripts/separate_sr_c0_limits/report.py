# -*- coding: utf-8 -*-
# Read-only native-limit report with physical C0 continuum continuity per era.
from pathlib import Path
import argparse,csv,datetime,hashlib,html,json,math,os
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_separate_sr_c0_limits'
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
    fig.suptitle('Continuità fisica C0 del continuo al bordo SR/CR, per era (C1 non imposta)\nQuattro ere combinate | limiti 95% CLs su r | μ Υ libere, sigma/code fisse',fontsize=13)
    provisional=any(q['six_quantiles'] and not q['expected_bands_verified'] for q in rows)
    note='Bande native provvisorie/non validate: valori conservati senza correzioni.' if provisional else 'Quantili nativi conservati; nessuna interpolazione delle masse mancanti.'
    fig.text(.5,.035,note+'\n× rosse: crossing osservato non verificato / errore; curva interrotta in quei punti.\nLa stessa CR compare nelle due analisi separate: queste curve non sono un limite combinato.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.81,bottom=.22,wspace=.24)
    for ext in ['png','pdf']:
        target=R/f'limits.{ext}';temp=R/f'.limits_{os.getpid()}.{ext}';fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(R/'limits_plot_cache.json',dict(signature=signature))


def comparison_plot(current_rows,points):
    previous_root=P/'runs/20261031_separate_sr_cheb_limits'
    if not any(finite(q['observed']) or finite(q['expected_50']) for q in current_rows):return
    current={(q['region'],q['mass']):q for q in current_rows}
    rows=[]
    for region,mass in points:
        path=previous_root/region/f'm{mass:.1f}/results.json';old=read(path,{})
        limits=old.get('limits') or {};receipt=old.get('receipt') or {};check=old.get('observed_limit_check') or {}
        new=current.get((region,mass),{})
        rows.append(dict(region=region,mass=mass,
            previous_observed=limits.get('observed'),previous_expected_50=limits.get('expected_50'),
            previous_observed_included=check.get('verified') is True and receipt.get('exit_code')==0,
            previous_flags=old.get('flags',[]),previous_source=str(path),
            C0_observed=new.get('observed'),C0_expected_50=new.get('expected_50'),
            C0_observed_included=new.get('observed_curve_included',False),C0_flags=new.get('flags',[]),
            C0_source=str(R/new['source']) if new else None))
    metadata=dict(previous_model='Separate SR+CR continua without an equality constraint at their shared physical boundary.',
        C0_model='Equality of physical continuum intensity at the SR/CR boundary per era. No equality of first derivatives is imposed.',
        shared_settings='Chebyshev6 low SR or Chebyshev4 high SR, Chebyshev4 CR, eight floating Upsilon means, fixed sigma/tails, four eras combined.',
        observed_policy='Only native observed values with verified crossing and exit code zero enter the curves. Other stored values remain in this data table but form gaps in the curves.',
        expected_policy='Native medians shown without bands for readability; expected bands/medians are not independently validated.',
        statistical_relation='Low/high analyses share CR data and are correlated; results are not combined or selected according to their observed limit.',rows=rows)
    signature=hashlib.sha256(json.dumps(clean(dict(metadata=metadata,plot_version=1)),sort_keys=True).encode()).hexdigest()
    cache=R/'comparison_c0_vs_unconstrained_cache.json'
    if read(cache,{}).get('signature')==signature and all((R/f'comparison_c0_vs_unconstrained.{ext}').exists() for ext in ['png','pdf','json','csv']):return
    write_json(R/'comparison_c0_vs_unconstrained.json',metadata);csv_write(R/'comparison_c0_vs_unconstrained.csv',rows)
    fig,axes=plt.subplots(1,2,figsize=(14,6.8))
    for region,ax in zip(REGIONS,axes):
        region_rows=[q for q in rows if q['region']==region]
        lookup={q['mass']:q for q in region_rows}
        requested=[q['mass'] for q in region_rows]
        if not requested:continue
        grid=[n/10 for n in range(round(min(requested)*10),round(max(requested)*10)+1)]
        for prefix,color,label in [('previous','#0072B2','Senza raccordo'),('C0','#7B3294','Con C0')]:
            observed=[lookup.get(m,{}).get(prefix+'_observed',np.nan) if lookup.get(m,{}).get(prefix+'_observed_included',False) else np.nan for m in grid]
            expected=[lookup.get(m,{}).get(prefix+'_expected_50',np.nan) for m in grid]
            ax.plot(grid,observed,'o-',ms=3,color=color,lw=1.6,label=label+': osservato verificato')
            ax.plot(grid,expected,'--',marker='_',ms=5,color=color,lw=1.4,label=label+': mediano atteso')
        ax.axvspan(9.093333333333334,10.506666666666668,color='#888888',alpha=.1,label='Ipotesi nella CR')
        n=sum(finite(q['C0_expected_50']) for q in region_rows)
        ax.set_title(f'{LABELS[region]} Chebyshev{DEGREES[region]} + CR Chebyshev4 | C0 {n}/{len(requested)}')
        ax.set_yscale('log');ax.set_xlim(min(grid)-.06,max(grid)+.06);ax.set_xlabel('Massa ipotizzata [GeV]')
        ax.set_ylabel('Limite superiore 95% CLs su r');ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best')
    fig.suptitle('Confronto di modelli: senza raccordo e con continuità fisica C0 al bordo SR/CR\nQuattro ere combinate; continuità delle derivate C1 non imposta',fontsize=13)
    fig.text(.5,.035,'Osservati: solo crossing verificati; le curve si interrompono ai valori mancanti o non verificati.\nMediani attesi nativi non validati; bande omesse per leggibilità. Le analisi bassa/alta sono separate e correlate tramite la CR.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.81,bottom=.2,wspace=.24)
    for ext in ['png','pdf']:
        target=R/f'comparison_c0_vs_unconstrained.{ext}';temp=R/f'.comparison_c0_{os.getpid()}.{ext}';fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(cache,dict(signature=signature))



def collect_fits(points):
    rows=[]
    for region,mass in points:
        path=R/region/f'm{mass:.1f}/fit_results.json';q=read(path)
        if q:
            q=dict(q,region=region,source=str(path.relative_to(R)));rows.append(q)
    return rows


def fit_gallery(q):
    region=q['region'];mass=q['mass'];directory=R/region/f'm{mass:.1f}'
    data_path=directory/'bins.json';result_path=directory/'fit_results.json'
    if not data_path.exists():return None
    target=R/'fits'/region/f'm{mass:.1f}';target.mkdir(parents=True,exist_ok=True)
    stamp=dict(version=2,bins_mtime=data_path.stat().st_mtime_ns,result_mtime=result_path.stat().st_mtime_ns)
    outputs={ext:str((target/f'fit_full_log.{ext}').relative_to(R)) for ext in ['png','pdf']}
    if read(target/'plot_cache.json')==stamp and all((R/value).exists() for value in outputs.values()):return outputs
    curves=read(data_path,{})
    if not all(f'{part}_ch{i}' in curves.get('s',{}) for part in ['SR','CR'] for i in range(1,5)):return None
    parts=['SR','CR'] if region=='low' else ['CR','SR']
    boundary=9.093333333333334 if region=='low' else 10.506666666666668
    domain=(6,10.506666666666668) if region=='low' else (9.093333333333334,12)
    fig=plt.figure(figsize=(14,10));outer=fig.add_gridspec(2,2,hspace=.32,wspace=.2)
    negative_signal=False;eras=['2022','2022EE','2023','2023BPix']
    for i,(era,cell) in enumerate(zip(eras,outer),1):
        sub=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.06)
        ax,rx=fig.add_subplot(sub[0]),fig.add_subplot(sub[1]);ymax=1
        for part in parts:
            channel=f'{part}_ch{i}';item=curves['s'][channel]
            x=np.asarray(item['mass'],dtype=float);data=np.asarray(item['data'],dtype=float)
            pred=np.asarray(item['total'],dtype=float);width=float(item['bin_width'])
            edges=np.r_[x-width/2,x[-1]+width/2];first=part==parts[0]
            ymax=max(ymax,float(np.max(data+np.sqrt(np.maximum(data,0)))))
            pos=data>0
            ax.errorbar(x[pos],data[pos],yerr=np.sqrt(data[pos]),fmt='.',color='black',ms=1.8,lw=.4,label='Dati' if first else None)
            ax.stairs(np.where(pred>0,pred,np.nan),edges,baseline=None,color='#009E73',lw=1.6,label='Totale S+B' if first else None)
            b=curves.get('b',{}).get(channel)
            if b:
                bv=np.asarray(b['total'],dtype=float)
                ax.stairs(np.where(bv>0,bv,np.nan),edges,baseline=None,color='gray',ls='--',lw=1,label='Solo fondo' if first else None)
            for name,raw in item['components'].items():
                values=np.asarray(raw,dtype=float);is_continuum=False
                if 'upsilon1s' in name:color,label='#E69F00','Υ(1S)'
                elif 'upsilon2s' in name:color,label='#CC79A7','Υ(2S)'
                elif 'Zd' in name:
                    color,label='#D55E00','Segnale DCB (parte positiva)';negative_signal|=bool(np.any(values<0))
                else:
                    color,label=('#0072B2','Continuo SR') if part=='SR' else ('#56B4E9','Continuo CR')
                    is_continuum=True
                ax.stairs(np.where(values>0,values,np.nan),edges,baseline=None,color=color,lw=1.1,label=label if first or is_continuum else None)
            residual=np.full_like(pred,np.nan);ok=np.isfinite(pred)&(pred>0)
            residual[ok]=(data[ok]-pred[ok])/np.sqrt(pred[ok]);rx.plot(x,residual,'.',color='#009E73',ms=2)
        for axis in [ax,rx]:
            axis.axvspan(9.093333333333334,10.506666666666668,color='#999999',alpha=.12)
            axis.axvline(boundary,color='#555555',ls=':',lw=.9)
            axis.set_xlim(*domain);axis.grid(alpha=.12)
        ax.set_title(era);ax.set_yscale('log');ax.set_ylim(.15,ymax*1.9)
        ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False)
        rx.axhline(0,color='gray',lw=.7);rx.set_ylabel('(D−M)/√M');rx.set_xlabel('Massa [GeV]')
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,ncol=4,fontsize=8.3,loc='upper center',bbox_to_anchor=(.5,.925),frameon=False)
    fig.suptitle(f'{LABELS[region]} + CR | mX={mass:.1f} GeV | Chebyshev{DEGREES[region]} / Chebyshev4\nContinuità C0 dell’intensità fisica al bordo {boundary:.6f} GeV; C1 non imposta; otto μ Υ libere',fontsize=13)
    stages=q.get('fit',{});b=stages.get('b',{});sb=stages.get('s',{})
    errors=[a.get('continuity',{}).get('maximum_relative_difference') for a in [b,sb]]
    checked=all(a.get('continuity_passed') is True for a in [b,sb])
    c0='PASS' if checked else 'FLAG / non disponibile'
    err=max(errors) if all(finite(value) for value in errors) else None
    invalid=q.get('point_valid') is not True
    flag_note='PUNTO S+B NON VALIDO' if invalid else 'Punto S+B valido; covarianza da verificare' if q.get('numerical_ok') is not True else 'Punto S+B valido'
    fig.text(.5,.075,f'{flag_note} | r̂={fmt(q.get("r"))}; B/S+B status={b.get("status","—")}/{sb.get("status","—")}, covQual={b.get("covQual","—")}/{sb.get("covQual","—")}, EDM={fmt(b.get("edm"))}/{fmt(sb.get("edm"))}',ha='center',fontsize=8.5,color='#C00000' if invalid else '#192639')
    fig.text(.5,.055,f'C0 post-fit {c0}, max scarto relativo={fmt(err)}; flag sorgente: '+('; '.join(q.get('flags',[])) or 'nessuno'),ha='center',fontsize=8.2,color='#C00000' if invalid else '#192639')
    note='Fascia grigia: CR. Curve native disegnate per regione; continuo SR e CR accoppiati nel valore al bordo.'
    lognote='Componenti non positive non visualizzate in log; tutti i bin restano nei residui.'
    if negative_signal:lognote+=' DCB stimata negativa.'
    fig.text(.5,.032,note,ha='center',fontsize=8.5);fig.text(.5,.012,lognote,ha='center',fontsize=8)
    fig.subplots_adjust(top=.85,bottom=.14)
    for ext,relative in outputs.items():
        out=R/relative;tmp=out.with_name(f'.fit_full_log_{os.getpid()}.{ext}');fig.savefig(tmp,dpi=145);tmp.replace(out)
    plt.close(fig);write_json(target/'plot_cache.json',stamp)
    return outputs


def fit_summary_row(q):
    row=dict(region=q['region'],mass=q['mass'],status=q.get('status'),r=q.get('r'),r_error=q.get('r_error'),
        point_valid=q.get('point_valid'),numerical_ok=q.get('numerical_ok'),flags=q.get('flags',[]),source=q['source'])
    for name in ['b','s']:
        stage=q.get('fit',{}).get(name,{})
        for key in ['status','covQual','edm','point_valid','numerical_ok','mean_boundary','continuity_passed']:row[f'{name}_{key}']=stage.get(key)
        row[f'{name}_C0_max_relative_difference']=stage.get('continuity',{}).get('maximum_relative_difference')
    return row



def collect_significance(points):
    rows=[]
    for region,mass in points:
        path=R/region/f'm{mass:.1f}/significance_results.json';q=read(path)
        if not q:continue
        receipt=q.get('receipt') or {};reference=q.get('reference_FD') or {};agreement=q.get('FD_agreement') or {}
        native_finite=finite(q.get('Z')) and finite(q.get('p0')) and 0<=q['p0']<=1
        complete=q.get('status')=='complete' and receipt.get('exit_code')==0 and native_finite
        rows.append(dict(region=region,mass=mass,Z=q.get('Z'),p0=q.get('p0'),log10_p0=q.get('log10_p0'),p0_underflow=q.get('p0_underflow'),status=q.get('status'),
            exit_code=receipt.get('exit_code'),native_finite=native_finite,method_complete=complete,
            numerical_ok=q.get('numerical_ok'),plot_included=complete and q.get('numerical_ok') is True,
            flags=q.get('flags',[]),covariance_ok=q.get('covariance_ok'),reference_FD=reference,FD_agreement=agreement,
            native_output=q.get('native_output'),source=str(path.relative_to(R))))
    return rows


def significance_plot(rows,points):
    if not any(q['native_finite'] for q in rows):return
    signature=hashlib.sha256(json.dumps(clean(dict(rows=rows,points=points,plot_version=3)),sort_keys=True).encode()).hexdigest()
    cache=R/'pvalues_plot_cache.json'
    if read(cache,{}).get('signature')==signature and all((R/f'pvalues.{ext}').exists() for ext in ['png','pdf']):return
    fig,axes=plt.subplots(1,2,figsize=(14,6.8))
    for region,ax,color in zip(REGIONS,axes,['#0072B2','#7B3294']):
        requested=[mass for reg,mass in points if reg==region];subset=[q for q in rows if q['region']==region]
        if not requested:continue
        lookup={q['mass']:q for q in subset};grid=[n/10 for n in range(round(min(requested)*10),round(max(requested)*10)+1)]
        vals=np.array([max(lookup[m]['p0'],1e-300) if m in lookup and lookup[m]['plot_included'] else np.nan for m in grid])
        ax.plot(grid,vals,'o-',ms=3,color=color,lw=1.5,label='p₀ locale nativo: controlli superati')
        bad=[q for q in subset if q['native_finite'] and not q['plot_included']]
        positive=[max(q['p0'],1e-300) for q in subset if q['plot_included']]
        bottom=max(1e-301,min([1e-4]+[p*.35 for p in positive]))
        visible_bad=[q for q in bad if q['p0']>=bottom]
        clipped_bad=[q for q in bad if q['p0']<bottom]
        if visible_bad:ax.plot([q['mass'] for q in visible_bad],[q['p0'] for q in visible_bad],'x',ms=8,mew=1.7,color='#C00000',label='Flag numerici / riferimento FD non valido',zorder=5)
        if clipped_bad:ax.plot([q['mass'] for q in clipped_bad],[bottom*1.2]*len(clipped_bad),'v',ms=8,color='#C00000',label=f'Flag: p₀ sotto scala (< {bottom:.1g})',zorder=5)
        for z in range(1,6):
            pz=.5*math.erfc(z/math.sqrt(2))
            if pz>=bottom:
                ax.axhline(pz,color='gray',ls=':',lw=.8);ax.text(.985,pz,f'{z}σ',ha='right',va='bottom',fontsize=8,transform=ax.get_yaxis_transform(),color='#555555')
        ax.axvspan(9.093333333333334,10.506666666666668,color='#888888',alpha=.1,label='Ipotesi nella CR')
        ax.set_title(f'{LABELS[region]} Chebyshev{DEGREES[region]} + CR Chebyshev4\n{sum(q["method_complete"] for q in subset)}/{len(requested)} esiti nativi')
        ax.set_yscale('log');ax.set_ylim(bottom,1.1);ax.set_xlim(min(grid)-.06,max(grid)+.06)
        ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('p-value locale p₀');ax.grid(alpha=.15);ax.legend(fontsize=8,loc='lower left',bbox_to_anchor=(0,.12))
    fig.suptitle('Significance nativa: continuità fisica C0 SR/CR, quattro ere combinate\nSignificanze locali asintotiche; nessuna correzione look-elsewhere',fontsize=13)
    note='× rosse: controllo numerico non superato; curva interrotta. ▼ rossi: valori segnalati sotto la scala del pannello.\nValori nativi integri in CSV/JSON. FD serve soltanto per controllo, mai per sostituire Z o p₀. Analisi correlate tramite la CR.'
    if any(q['p0']==0 for q in rows if q['native_finite']):note+='\np₀ nativo nullo per underflow è conservato come 0 nei dati; i punti segnalati restano marcati sotto scala.'
    fig.text(.5,.035,note,ha='center',fontsize=8.6);fig.subplots_adjust(top=.78,bottom=.21,wspace=.24)
    for ext in ['png','pdf']:
        target=R/f'pvalues.{ext}';temp=R/f'.pvalues_{os.getpid()}.{ext}';fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(cache,dict(signature=signature))


def report(manifest_path):
    R.mkdir(parents=True,exist_ok=True)
    manifest,points=load_manifest(manifest_path);rows=collect();lookup={(q['region'],q['mass']):q for q in rows}
    selected=[q for q in rows if (q['region'],q['mass']) in points]
    extras=[q for q in rows if (q['region'],q['mass']) not in points]
    plot(selected,points);comparison_plot(selected,points)
    significance_requested=bool(manifest and (manifest.get('significance_requested') or 'Significance' in manifest.get('requested_methods',[])))
    significance_rows=collect_significance(points) if significance_requested else []
    significance_lookup={(q['region'],q['mass']):q for q in significance_rows}
    significance_plot(significance_rows,points)
    fits_requested=bool(manifest and manifest.get('fits_requested'))
    fit_records=collect_fits(points) if fits_requested else []
    fit_lookup={(q['region'],q['mass']):q for q in fit_records}
    fit_rows=[];render_errors=[]
    priority={('low',7.5):0,('high',11.):1}
    ordered=sorted(fit_records,key=lambda q:(priority.get((q['region'],q['mass']),2),REGIONS.index(q['region']),q['mass']))
    for i,q in enumerate(ordered,1):
        row=fit_summary_row(q)
        try:
            row['gallery']=fit_gallery(q)
            if not row['gallery']:render_errors.append(dict(region=q['region'],mass=q['mass'],error='Missing or incomplete binned predictions'))
        except Exception as exc:
            row['gallery']=None;render_errors.append(dict(region=q['region'],mass=q['mass'],error=repr(exc)))
        fit_rows.append(row)
        if i<=2 or i%10==0 or i==len(ordered):
            print(json.dumps(dict(rendered=i,target=len(ordered),last_region=q['region'],last_mass=q['mass'],render_errors=len(render_errors))),flush=True)
    fit_rows.sort(key=lambda q:(REGIONS.index(q['region']),q['mass']))
    fit_plot_lookup={(q['region'],q['mass']):q for q in fit_rows}
    counts={region:dict(target=sum(reg==region for reg,m in points),results=sum(q['region']==region for q in selected),limits=sum(q['region']==region and q['six_quantiles'] for q in selected),verified_crossings=sum(q['region']==region and q['observed_curve_included'] for q in selected),completed_methods=sum(q['region']==region and q['method_complete'] for q in selected),unvalidated_observed=sum(q['region']==region and finite(q['observed']) and not q['observed_curve_included'] for q in selected)) for region in REGIONS}
    for region in REGIONS:
        subset=[q for q in fit_rows if q['region']==region]
        counts[region].update(fit_results=len(subset),fit_methods_complete=sum(q['status']=='complete' for q in subset),fit_galleries=sum(bool(q['gallery']) for q in subset),fit_points_valid=sum(q['point_valid'] is True for q in subset),fit_points_invalid=sum(q['point_valid'] is not True for q in subset),fit_covariance_ok=sum(q['numerical_ok'] is True for q in subset),fit_s_covQual2=sum(q['s_covQual']==2 for q in subset),fit_b_valid=sum(q['b_point_valid'] is True for q in subset),postfit_C0_passed=sum(q['b_continuity_passed'] is True and q['s_continuity_passed'] is True for q in subset))
    all_limits=bool(manifest) and not manifest.get('CR_masses_pending',False) and all(lookup.get(point,{}).get('method_complete',False) for point in points)
    all_fits=not fits_requested or all(fit_lookup.get(point,{}).get('status')=='complete' and fit_plot_lookup.get(point,{}).get('gallery') for point in points)
    all_significance=not significance_requested or all(significance_lookup.get(point,{}).get('method_complete',False) for point in points)
    completed=bool(all_limits and all_fits and all_significance and not render_errors)
    limit_outcomes=sum(any((R/reg/f'm{m:.1f}'/name).exists() for name in ['completion.json','worker_failure.json']) for reg,m in points)
    fit_outcomes=sum(any((R/reg/f'm{m:.1f}'/name).exists() for name in ['fits_completion.json','fit_worker_failure.json']) for reg,m in points)
    remaining=read(R/'remaining_runner.json',{});fit_runner=read(R/'fit_runner.json',{})
    production_closed=remaining.get('status')=='completed' and limit_outcomes==len(points) and (not fits_requested or (fit_runner.get('status')=='completed' and fit_outcomes==len(points)))
    significance_outcomes=sum(q.get('status') in ['complete','failed'] for q in significance_rows)
    production_closed=production_closed and (not significance_requested or significance_outcomes==len(points))
    for region in REGIONS:
        sr=[q for q in significance_rows if q['region']==region]
        counts[region].update(significance_outcomes=len(sr),significance_complete=sum(q['method_complete'] for q in sr),significance_checked=sum(q['plot_included'] for q in sr),significance_flagged=sum(q['native_finite'] and not q['plot_included'] for q in sr))
    stopped=(R/'stop_record.json').exists()
    watcher_done=completed or stopped or (production_closed and not render_errors)
    now=datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds')
    out=dict(updated_at=now,completed=completed,production_closed=production_closed,stopped_by_user=stopped,watcher_done=watcher_done,renderer_parent_pid=os.getppid(),manifest_available=manifest is not None,manifest_path=str(manifest_path),target=len(points),counts=counts,points=[dict(region=r,mass=m) for r,m in points],results=selected,unrequested_results=extras,fits_requested=fits_requested,fit_outcomes=fit_outcomes,limit_outcomes=limit_outcomes,render_errors=render_errors,significance_requested=significance_requested,significance_outcomes=significance_outcomes,
        interpretation='Physical C0 continuity of the continuum intensity is imposed at each SR/CR boundary per era; first-derivative C1 continuity is not imposed. Separate low-SR+CR and high-SR+CR analyses, statistically correlated through their shared CR data. The common CR is not multiplied or combined. Overlapping CR mass hypotheses are shown for both declared likelihoods, without selecting a curve based on its observed limit. Native limits and bands are unchanged; missing masses are not interpolated. Completion accounts for method outcomes and figures, and does not certify fit validity or expected-band validation.')
    write_json(R/'limits.json',out);csv_write(R/'limits.csv',selected)
    write_json(R/'fits.json',dict(updated_at=now,requested=fits_requested,target=len(points),results=fit_rows,render_errors=render_errors,interpretation='Native B and S+B predictions. Point validity, covariance quality, exact post-fit continuity and method completion are distinct. Residuals (data-model)/sqrt(model) are visual diagnostics, not an additional calibrated goodness-of-fit test.'))
    csv_write(R/'fits.csv',fit_rows)
    if significance_requested:
        write_json(R/'pvalues.json',dict(updated_at=now,target=len(points),results=significance_rows,interpretation='Native Combine Significance Z and p0 are unchanged. Local asymptotic significance, no look-elsewhere correction. FitDiagnostics is only a numerical cross-check. Correlated low/high analyses share CR data.'))
        csv_write(R/'pvalues.csv',significance_rows)
    write_json(R/'progress.json',{k:out[k] for k in ['updated_at','completed','production_closed','stopped_by_user','watcher_done','renderer_parent_pid','manifest_available','target','counts','fits_requested','fit_outcomes','limit_outcomes','render_errors','significance_requested','significance_outcomes']})
    tables=[];galleries=[]
    for region in REGIONS:
        cells=[];cards=[]
        for reg,m in points:
            if reg!=region:continue
            q=lookup.get((reg,m),{})
            flags=html.escape('; '.join(q.get('flags',[]))) or '—'
            source=f'<a href="{q["source"]}">JSON</a>' if q else 'Non disponibile' if production_closed or stopped else 'In attesa'
            crossing='Verificato' if q.get('observed_curve_included') else 'NON VERIFICATO: valore non validato' if q else '—'
            row_style=' class="unvalidated"' if q and not q.get('observed_curve_included') else ''
            method_label=html.escape(str(q.get('status') or '—'))+'; exit '+str(q.get('exit_code','—'))
            cells.append(f'<tr{row_style}><td>{m:.1f}</td><td>{fmt(q.get("observed"))}</td><td>{fmt(q.get("expected_2p5"))}</td><td>{fmt(q.get("expected_16"))}</td><td>{fmt(q.get("expected_50"))}</td><td>{fmt(q.get("expected_84"))}</td><td>{fmt(q.get("expected_97p5"))}</td><td>{crossing}</td><td>{method_label}</td><td>{flags}</td><td>{source}</td></tr>')
            if not fits_requested:continue
            f=fit_plot_lookup.get((reg,m),{});gallery=f.get('gallery');original=fit_lookup.get((reg,m),{})
            fs=original.get('fit',{});b=fs.get('b',{});s=fs.get('s',{})
            flag_text=html.escape('; '.join(f.get('flags',[]))) or 'Nessuno'
            badge='Punto S+B valido' if f.get('point_valid') is True else 'PUNTO S+B NON VALIDO' if f else 'Fit non disponibile'
            covariance='Covarianza B/S+B di qualità 3' if f.get('numerical_ok') is True else 'Controllare covarianza / flag: distinta dalla validità del punto'
            style=' invalid' if f and f.get('point_valid') is not True else ''
            thumbnail=f'<a href="{gallery["png"]}"><img loading="lazy" src="{gallery["png"]}" alt="Fit {LABELS[reg]} massa {m:.1f}, quattro ere"></a><p><a href="{gallery["png"]}">PNG</a> · <a href="{gallery["pdf"]}">PDF</a> · <a href="{f["source"]}">Fit originale JSON</a> · <a href="{reg}/m{m:.1f}/bins.json">Bin e componenti</a></p>' if gallery else '<p>Figura non disponibile: vedere esiti del rendering.</p>'
            stage_rows=''.join(f'<tr><td>{label}</td><td>{stage.get("status","—")}</td><td>{stage.get("covQual","—")}</td><td>{fmt(stage.get("edm"))}</td><td>{"valido" if stage.get("point_valid") is True else "NON VALIDO"}</td><td>{"PASS" if stage.get("continuity_passed") is True else "FLAG"}; {fmt(stage.get("continuity",{}).get("maximum_relative_difference"))}</td></tr>' for label,stage in [('B',b),('S+B',s)])
            boundaries=html.escape(json.dumps({stage:fs.get(stage,{}).get('mean_boundary',[]) for stage in ['b','s']}))
            cards.append(f'<article class="fitcard{style}" id="fit-{reg}-{m:.1f}"><h3>{LABELS[reg]}, m = {m:.1f} GeV</h3><p><b>{badge}</b> · r̂ = {fmt(f.get("r"))} ± {fmt(f.get("r_error"))}<br>{covariance}</p>{thumbnail}<div class="scroll"><table><thead><tr><th>Fit</th><th>Status</th><th>covQual</th><th>EDM</th><th>Punto</th><th>C0: max scarto</th></tr></thead><tbody>{stage_rows}</tbody></table></div><p><b>Flag sorgente:</b> {flag_text}<br>μ ai bordi dei range: <code>{boundaries}</code></p></article>')
        tables.append(f'<section><h2>Valori dei limiti: {LABELS[region]} + CR</h2><div class="scroll"><table><thead><tr><th>Massa</th><th>Osservato</th><th>Atteso −2σ</th><th>−1σ</th><th>Mediano</th><th>+1σ</th><th>+2σ</th><th>Crossing obs.</th><th>Esito metodo</th><th>Flag sorgente</th><th>Risultato</th></tr></thead><tbody>'+''.join(cells)+'</tbody></table></div></section>')
        if fits_requested:galleries.append(f'<section id="gallery-{region}"><h2>Fit: {LABELS[region]} + CR</h2><div class="gallery">'+''.join(cards)+'</div></section>')
    refresh='' if watcher_done else '<meta http-equiv="refresh" content="60">'
    plot_html='<p class="pending">In attesa dei primi limiti nativi.</p>'
    if (R/'limits.png').exists():plot_html='<a href="limits.png"><img src="limits.png" alt="Limiti separati SR bassa e alta con bande"></a><p><a href="limits.pdf">PDF del grafico</a></p>'
    comparison_html='<p class="pending">Il confronto apparirà ai primi risultati C0.</p>'
    if (R/'comparison_c0_vs_unconstrained.png').exists():comparison_html='<a href="comparison_c0_vs_unconstrained.png"><img src="comparison_c0_vs_unconstrained.png" alt="Confronto senza raccordo e con C0"></a><p><a href="comparison_c0_vs_unconstrained.pdf">PDF</a> · <a href="comparison_c0_vs_unconstrained.csv">Valori CSV</a> · <a href="comparison_c0_vs_unconstrained.json">Valori e provenienza JSON</a></p>'
    manifest_note='Elenco dei punti definito dal manifest.' if manifest else 'Manifest ancora in attesa; nessun completamento dichiarato prima del manifest.'
    links='<a href="limits.json">Dati e configurazione JSON</a> · <a href="progress.json">Avanzamento</a>'
    if selected:links+=' · <a href="limits.csv">Limiti CSV</a>'
    if fit_rows:links+=' · <a href="fits.csv">Riepilogo fit CSV</a> · <a href="fits.json">Fit e gallerie JSON</a>'
    for name,label in [('completion_audit.json','Audit indipendente di completezza fit e limiti'),('completion_audit.md','Audit fit e limiti leggibile'),('significance_completion.json','Audit significanze native'),('pilot_model_audit.json','Audit costruttivo piloti'),('all_mass_model_audit.json','Audit costruttivo di tutti i modelli'),('workspace_manifest.json','Provenienza e hash dei workspace')]:
        if (R/name).exists():links+=f' · <a href="{name}">{label}</a>'
    if manifest_path.exists():links+=f' · <a href="{os.path.relpath(manifest_path,R)}">Manifest</a>'
    significance_html=''
    if significance_requested:
        sigpic='<p class="pending">In attesa dei risultati nativi Significance.</p>'
        if (R/'pvalues.png').exists():sigpic='<a href="pvalues.png"><img src="pvalues.png" alt="P-value locali nativi con C0"></a><p><a href="pvalues.pdf">PDF</a> · <a href="pvalues.csv">CSV</a> · <a href="pvalues.json">JSON e provenienza</a></p>'
        sigcells=[]
        for reg,m in points:
            q=significance_lookup.get((reg,m),{});ref=q.get('reference_FD',{});agreement=q.get('FD_agreement',{})
            state='Controlli superati' if q.get('plot_included') else 'FLAG / non validato' if q else 'In attesa'
            source=f'<a href="{q["source"]}">JSON nativo</a>' if q else '—'
            sigflags=html.escape('; '.join(q.get('flags',[]))) or '—'
            sigcells.append(f'<tr class="{"unvalidated" if q and not q.get("plot_included") else ""}"><td>{LABELS[reg]}</td><td>{m:.1f}</td><td>{fmt(q.get("Z"))}</td><td>{fmt(q.get("p0"))}</td><td>{state}</td><td>{html.escape(str(ref.get("point_valid","—")))}</td><td>{fmt(agreement.get("absolute_Z_difference"))}</td><td>{sigflags}</td><td>{source}</td></tr>')
        significance_html=f'<section id="pvalues"><h2>P-value e significanza locale</h2><p>{sum(q["method_complete"] for q in significance_rows)}/{len(points)} esiti nativi; {sum(q["plot_included"] for q in significance_rows)} con controlli numerici superati. Z e p₀ provengono dal metodo Combine Significance: non sono ricavati dalle figure né sostituiti dal confronto con FitDiagnostics. I valori con controlli falliti restano visibili come croci rosse; i triangoli rossi indicano quelli sotto la scala del pannello. La tabella conserva i numeri nativi, incluso p₀ = 0 da underflow. Le significanze sono locali e asintotiche, senza correzione look-elsewhere; le due analisi sono correlate attraverso la CR.</p>{sigpic}<div class="scroll"><table><thead><tr><th>Analisi</th><th>Massa</th><th>Z nativo</th><th>p₀ nativo</th><th>Controlli</th><th>Riferimento FD valido</th><th>|ΔZ| FD</th><th>Flag originali</th><th>Sorgente</th></tr></thead><tbody>'+''.join(sigcells)+'</tbody></table></div></section>'
    fit_total=sum(c['fit_results'] for c in counts.values());valid_total=sum(c['fit_points_valid'] for c in counts.values());cov2=sum(c['fit_s_covQual2'] for c in counts.values());c0_total=sum(c['postfit_C0_passed'] for c in counts.values())
    fit_intro=f'<section><h2>Fit e componenti nelle quattro ere</h2><p><b>{fit_total}/{len(points)} esiti FitDiagnostics; {valid_total}/{fit_total} punti S+B validi; {cov2} covarianze S+B con covQual = 2.</b> La qualità della matrice di covarianza è riportata separatamente dalla validità della stima centrale. Continuità C0 esatta post-fit: {c0_total}/{fit_total} passano in entrambi i fit B e S+B.</p><p>Un grafico per ipotesi: quattro ere, dati, totale S+B, totale solo fondo e componenti S+B; scala logaritmica con residui (D−M)/√M. La fascia grigia indica la CR. I residui sono diagnostici visivi, non un nuovo test di goodness of fit. Le componenti non positive sono omesse soltanto dalla scala logaritmica. Le curve sono disegnate nei propri intervalli; il controllo C0 usa la densità esatta al bordo, non i centri di due bin adiacenti. Gli audit costruttivi e le verifiche post-fit sono distinti.</p><p><a href="#gallery-low">Galleria SR bassa</a> · <a href="#gallery-high">Galleria SR alta</a> · <a href="fits.csv">Tabella dei fit CSV</a> · <a href="fits.json">Esiti e link JSON</a></p></section>' if fits_requested else ''
    render_notice=f'<p class="notice">Errori rendering: {html.escape(json.dumps(render_errors))}</p>' if render_errors else ''
    status_text='Tutti gli esiti di calcolo e i grafici previsti sono disponibili; completamento distinto dalla validità dei fit e dei limiti.' if completed else 'Produzione interrotta: sono mostrati soltanto i risultati disponibili.' if stopped else 'Esecuzioni concluse con risultati mancanti o errori: vedere le tabelle.' if production_closed else 'Risultati in corso di produzione.'
    page=f'''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">{refresh}<title>Xee — fit, limiti e p-value con continuità C0 al bordo SR/CR</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1450px;margin:auto;padding:25px;background:#f4f6fa;color:#192639}}section,.card{{background:white;border:1px solid #d8e1ec;border-radius:12px;padding:20px;margin:16px 0}}.cards,.gallery{{display:grid;grid-template-columns:1fr 1fr;gap:15px}}.fitcard{{border:1px solid #d8e1ec;padding:12px;border-radius:8px;min-width:0;font-size:13px}}.invalid{{border:2px solid #bc2929;background:#fff4f1}}.count{{font-size:28px;font-weight:700}}.pending{{color:#68798d}}.notice{{background:#fff6e7;padding:12px;border-left:4px solid #c58724}}img{{max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:8px;text-align:left;border-bottom:1px solid #dde4ec}}th{{background:#edf2f8}}.scroll{{overflow:auto}}.unvalidated{{background:#fff0ec}}a{{color:#0758a2}}code{{overflow-wrap:anywhere}}@media(max-width:900px){{.cards,.gallery{{grid-template-columns:1fr}}body{{padding:12px}}}}</style></head><body>
<h1>Fit, limiti e p-value con continuità C0: SR bassa + CR e SR alta + CR</h1><p>{status_text}<br>Aggiornato: {now}</p>{render_notice}
<div class="cards"><div class="card"><div class="count">{counts['low']['limits']} / {counts['low']['target']}</div>SR bassa: Chebyshev6 + CR Chebyshev4<br>{counts['low']['verified_crossings']} crossing osservati verificati · {counts['low']['unvalidated_observed']} valori non validati<br>{counts['low']['fit_galleries']} figure fit · {counts['low']['fit_points_valid']} punti S+B validi</div><div class="card"><div class="count">{counts['high']['limits']} / {counts['high']['target']}</div>SR alta: Chebyshev4 + CR Chebyshev4<br>{counts['high']['verified_crossings']} crossing osservati verificati · {counts['high']['unvalidated_observed']} valori non validati<br>{counts['high']['fit_galleries']} figure fit · {counts['high']['fit_points_valid']} punti S+B validi</div></div>
<section><h2>Configurazione</h2><p>SR bassa: 6–9,093333 GeV. SR alta: 10,506667–12 GeV. CR: 9,093333–10,506667 GeV. Ogni curva deriva da una propria likelihood SR+CR: <b>la CR non viene moltiplicata e le due curve non sono combinate</b>. Le due analisi usano gli stessi dati della CR e sono statisticamente correlate.</p><p>Singole funzioni Chebyshev, grado 6 nella SR bassa e grado 4 nella SR alta; grado 4 nella CR in entrambe le analisi. <b>In ogni era, l’intensità fisica del continuo è uguale ai due lati del bordo SR/CR</b>: C0 a 9,093333 GeV per la SR bassa e a 10,506667 GeV per la SR alta. Nella parametrizzazione, N_CR è libero e N_SR = N_CR × (P_CR / I_CR) × (I_SR / P_SR), dove P è il valore al bordo e I l’integrale nel rispettivo intervallo. Si rimuovono quattro parametri liberi. Non è una prior gaussiana e non è imposta l’uguaglianza delle derivate (C1). Otto medie Υ libere, sigma e code fissate; normalizzazioni Υ condivise tra SR e CR; segnale DCB in entrambe le regioni selezionate con lo stesso r. {manifest_note} Le ipotesi 9,1–10,5 GeV nella CR sono valutate in entrambe le analisi; nessuna curva è scelta in base al limite osservato più favorevole.</p><p class="notice">Quattro ere combinate in ciascuna likelihood. Quantili e bande nativi conservati; gli eventuali flag di ricerca del bracket restano visibili nelle tabelle. Le bande senza verifica documentata restano provvisorie/non validate. Il crossing verificato non certifica l’intero minimo della likelihood. Sopra 11 GeV si usano le DCB parametriche e le normalizzazioni di segnale estrapolate nello studio già autorizzato, senza aggiungere una nuova incertezza di estrapolazione. Il completamento dei calcoli non equivale a validazione del modello.</p></section>
<section><h2>Limiti 95% CLs su r</h2><p>Osservato, atteso mediano e bande attese a 1 e 2 sigma. Le croci rosse riportano i valori nativi con crossing osservato non verificato o errore di esecuzione: non sono limiti validati. La curva osservata si interrompe in quei punti, mentre CSV e tabelle conservano integralmente numeri e flag originali. I punti mancanti restano mancanti. I pannelli usano la stessa unità r, con scale verticali logaritmiche proprie.</p>{plot_html}</section>
<section><h2>Confronto con il modello precedente senza raccordo</h2><p>Modelli distinti: precedente con continui SR/CR senza uguaglianza imposta al bordo; nuovo con continuità fisica C0 per era. Osservati tracciati soltanto dove il crossing è verificato, con interruzioni altrove. Mediani attesi nativi; bande omesse nel confronto per leggibilità. Le due analisi bassa/alta restano separate e correlate tramite i dati della CR.</p>{comparison_html}</section>
{significance_html}{fit_intro}{''.join(galleries)}{''.join(tables)}<section><h2>Dati e audit</h2><p>{links}</p><p>Risultati presenti fuori dall’elenco corrente: {len(extras)}. Il report legge i file sorgente senza modificarli.</p></section></body></html>'''
    temp=R/'index.html.tmp';temp.write_text(page,encoding='utf-8');temp.replace(R/'index.html')
    print(json.dumps(dict(updated_at=now,completed=completed,manifest_available=manifest is not None,target=len(points),counts=counts,render_errors=render_errors)),flush=True)
    return completed

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,default=R/'manifest.json');args=parser.parse_args();report(args.manifest)

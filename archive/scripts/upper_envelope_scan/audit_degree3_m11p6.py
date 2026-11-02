from pathlib import Path
import ROOT,json,math
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan/degree3';out={};keep=[]

def opened(mass):
    f=ROOT.TFile.Open(str(R/f'm{mass}/joint/workspace_seeded.root'));keep.append(f);return f.Get('w')
def rows(w):
    data=w.data('data_obs');a=[]
    for j in range(data.numEntries()):
        q=data.get(j);a.append((q.find('CMS_channel').getLabel(),q.getRealValue('mass'),q.getRealValue('mass_CR'),data.weight()))
    return sorted(a)
def bpred(w,mass):
    f=ROOT.TFile.Open(str(R/f'm{mass}/joint/fitdiagnostics/fitDiagnostics_upper.root'));keep.append(f);fit=f.Get('fit_b')
    for v in fit.floatParsFinal():w.var(v.GetName()).setVal(v.getVal())
    for era in ['2022','2022EE','2023','2023BPix']:
        bp=fit.floatParsFinal().find('local_b1_'+era);cp=fit.floatParsFinal().find('local_c1_'+era);idx=1 if bp and bp.getError()>0 and (not cp or cp.getError()==0) else 0;w.cat('pdfindex_'+era).setIndex(idx)
    w.var('r').setVal(0);pred=[]
    for reg,obsname in [('SR','mass'),('CR','mass_CR')]:
        for i in range(1,5):
            cp=w.obj('ModelConfig_bonly').GetPdf().getPdf(f'{reg}_ch{i}');pdf=next(q for q in cp.getComponents() if q.InheritsFrom('RooAddPdf'));x=w.var(obsname);obs=ROOT.RooArgSet(x);n=pdf.expectedEvents(obs);bw=(x.getMax()-x.getMin())/x.getBins()
            for j in range(x.getBins()):x.setVal(x.getMin()+(j+.5)*bw);pred.append(n*pdf.getVal(obs)*bw)
    return pred,dict(status=fit.status(),covQual=fit.covQual())
w=opened('11.6');other=opened('11.5');a=rows(w);b=rows(other);assert a==b
out.update(dataset_identical_to_mass11p5=True,entries=len(a),sum_weights=sum(q[3] for q in a))
out['channels']={}
for name in ['ModelConfig','ModelConfig_bonly']:
    pdf=w.obj(name).GetPdf();channels=[f'{r}_ch{i}' for r in ['SR','CR'] for i in range(1,5)];assert all(pdf.getPdf(ch) for ch in channels);out['channels'][name]=channels
out['signal_checks']=[];w.var('r').setVal(1)
for reg,sub,suffix in [('SR','signal',''),('CR','control','_CR')]:
    audit=json.loads((R/'m11.6'/sub/'audit.json').read_text())
    for era in audit['eras']:
        ref=era['processes']['Zd'];f=w.function('local_acceptance_Zd_'+era['era']+suffix);n=w.function('local_rate_Zd_'+era['era']+suffix)
        assert abs(f.getVal()-ref['fraction'])<1e-10
        assert abs(n.getVal()-ref['local_rate_r1'])<1e-7
        out['signal_checks'].append(dict(region=reg,era=era['era'],fraction=f.getVal(),rate=n.getVal(),reference=ref['local_rate_r1']))
p,q=bpred(w,'11.6');r,s=bpred(other,'11.5');maximum=max(abs(x-y)/max(1,abs(x)) for x,y in zip(p,r));out.update(B_quality_m11p6=q,B_quality_m11p5=s,max_relative_B_prediction_difference=maximum);assert maximum<1e-4 and q['status']==0 and q['covQual']==3
(R/'m11.6/preparation_semantic_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

"""Audit saved generation strengths and expected regional signal contributions."""
from pathlib import Path
import json,math
import ROOT
from generator_toys import R
from toys_helpers import atomic_json,channel_layout


def main():
    paths=list((R/'bias/generation').glob('m*/rExpected/truth.json'))
    paths += [p for p in [R/'bias/generation/common_r0/truth.json',R/'gof/generation/truth.json'] if p.exists()]
    rows=[];issues=[]
    for path in paths:
        truth=json.loads(path.read_text());file=ROOT.TFile.Open(str(path.parent/'generator.root'));w=file.Get('w')
        rtrue=float(truth['rtrue']);assert abs(w.var('r').getVal()-rtrue)<1e-12*max(1,abs(rtrue))
        channels=[]
        for channel,obsname in channel_layout(w):
            cp=w.obj('ModelConfig').GetPdf().getPdf(channel);pdf=next(p for p in cp.getComponents() if p.InheritsFrom('RooAddPdf'))
            x=w.var(obsname);obs=ROOT.RooArgSet(x);w.var('r').setVal(rtrue);expected=pdf.expectedEvents(obs)
            w.var('r').setVal(0);background=pdf.expectedEvents(obs);w.var('r').setVal(1);signal_r1=pdf.expectedEvents(obs)-background
            residual=expected-background-rtrue*signal_r1
            tolerance=1e-8*max(1,abs(expected),abs(background))
            if abs(residual)>tolerance:issues.append(dict(path=str(path),channel=channel,issue='injection_not_linear',residual=residual))
            mean=sum(a['counts'][channel] for a in truth['toy_audit'])/truth['ntoys']
            mc_se=math.sqrt(expected/truth['ntoys']) if expected>0 else None
            channels.append(dict(channel=channel,expected=expected,background=background,signal_r1=signal_r1,
                                 injected_signal=rtrue*signal_r1,injection_residual=residual,toy_mean=mean,
                                 toy_mean_z=(mean-expected)/mc_se if mc_se else None))
        file.Close();rows.append(dict(path=str(path),rtrue=rtrue,toys=truth['ntoys'],channels=channels))
    out=dict(ensembles=len(rows),rows=rows,issues=issues,passed=not issues,
             meaning='Saved r matches requested truth; all regional signal yields scale with that common r. Poisson sample means shown diagnostically, without selection.')
    atomic_json(R/'generation_strength_audit.json',out)
    print(json.dumps(dict(ensembles=len(rows),passed=out['passed'],max_abs_toy_mean_z=max(abs(c['toy_mean_z']) for row in rows for c in row['channels'] if c['toy_mean_z'] is not None)),indent=2))


if __name__=='__main__':main()

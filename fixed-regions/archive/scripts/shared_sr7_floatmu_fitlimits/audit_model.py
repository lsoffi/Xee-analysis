"""Read-only model-delta audit; no fit or minimization is performed."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import struct
import ROOT
import numpy as np
from model import P,R,ERAS,REGIONS


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(mass,deep=False):
    mass=f'{float(mass):.1f}';path=R/f'observed/m{mass}/joint/workspace_seeded.root'
    before=sha(path);handle=ROOT.TFile.Open(str(path));w=handle.Get('w')
    meta=json.loads(path.with_suffix('.audit.json').read_text())
    source=ROOT.TFile.Open(meta['source']);old=source.Get('w')
    assert sha(meta['source'])==meta['source_sha256']
    assert meta['degree']==7 and not meta['means_fixed']
    means=[v for v in w.allVars() if v.GetName().startswith('Ymean_')]
    assert len(means)==8 and all(not v.isConstant() for v in means)
    assert w.genobj('discreteParams').getSize()==4
    multi=[p for p in w.obj('ModelConfig').GetPdf().getComponents() if p.InheritsFrom('RooMultiPdf')]
    assert len(multi)==8 and abs(2*sum(p.getCorrection()for p in multi)-32)<1e-12
    assert w.data('data_obs').numEntries()==1800 and w.data('data_obs').sumEntries()==old.data('data_obs').sumEntries()
    assert w.obj('ModelConfig').GetGlobalObservables().getSize()==10
    assert w.obj('ModelConfig').GetPdf()._extraConstraints.getSize()==18
    rate_partitions=[];intensity_errors=[];mu_checks=[]
    for i,era in enumerate(ERAS,1):
        for proc in ['upsilon1s','upsilon2s']:
            mv=w.var(f'Ymean_{proc}_{era}');original=old.pdf(f'shapeBkg_{proc}_ch{i}')
            assert abs(mv.getMin()-(original.x0().getVal()-.2))<1e-12
            assert abs(mv.getMax()-(original.x0().getVal()+.2))<1e-12
            assert w.function(f'Ysigma_{proc}_{era}').isConstant()
            assert abs(w.function(f'Ysigma_{proc}_{era}').getVal()-original.sigmaL().getVal())<1e-12
            assert not w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}').isConstant()
            for j,p in enumerate([original.alphaL(),original.nL(),original.alphaR(),original.nR()]):
                assert w.function(f'Ytail_{proc}_{era}_{j}').isConstant()
                assert w.function(f'Ytail_{proc}_{era}_{j}').getVal()==p.getVal()
            assert all(w.pdf(f'shape_{proc}_{era}_{reg}').dependsOn(mv)for reg in REGIONS)
        idx=w.cat('pdfindex_'+era);saved=idx.getIndex()
        for family in [0,1]:
            idx.setIndex(family)
            total=w.var('local_yield_'+era).getVal()
            rates={reg:w.function(f'local_rate_dy_{era}_{reg}').getVal()for reg in ['SRL','SRH']}
            err=abs(sum(rates.values())/total-1);rate_partitions.append(err);assert err<1e-12
            if deep:
                c=[1.]+[w.var(f'local_c{j}_{era}').getVal()for j in range(1,8)]
                b=[1.]+[w.var(f'local_b{j}_{era}').getVal()for j in range(1,8)]
                integral=sum(w.function(f'integral_selected_{era}_{reg}').getVal()for reg in ['SRL','SRH'])
                for reg in ['SRL','SRH']:
                    x=w.var('mass_'+reg);norm=ROOT.RooArgSet(x);pdf=w.pdf(f'local_envelope_{era}_{reg}')
                    lo,hi,_=REGIONS[reg]
                    for frac in [.17,.51,.83]:
                        xx=lo+frac*(hi-lo);x.setVal(xx)
                        raw=np.polynomial.chebyshev.chebval((xx-9)/3,c) if family==0 else sum(b[j]*math.comb(7,j)*((xx-6)/6)**j*((12-xx)/6)**(7-j)for j in range(8))
                        expected=total*raw/integral;native=rates[reg]*pdf.getVal(norm)
                        err=abs(native/expected-1);intensity_errors.append(err);assert err<1e-9
        idx.setIndex(saved)
        for proc in ['Zd','upsilon1s','upsilon2s']:
            fractions=[w.function(f'local_acceptance_{proc}_{era}_{reg}').getVal()for reg in REGIONS]
            assert abs(sum(fractions)-1)<1e-10
        if deep:
            for proc in ['upsilon1s','upsilon2s']:
                mv=w.var(f'Ymean_{proc}_{era}');central=mv.getVal()
                x=ROOT.RooRealVar(f'audit_mass_{proc}_{era}','',9.,6.,12.)
                sig=w.function(f'Ysigma_{proc}_{era}')
                tails=[w.function(f'Ytail_{proc}_{era}_{j}')for j in range(4)]
                full=ROOT.RooCrystalBall(f'audit_full_{proc}_{era}','',x,mv,sig,sig,*tails)
                original_fractions={reg:w.function(f'local_acceptance_{proc}_{era}_{reg}').getVal()for reg in REGIONS}
                for delta in [-.03,.03]:
                    mv.setVal(min(mv.getMax()-.001,max(mv.getMin()+.001,central+delta)))
                    fractions=[]
                    for reg,(lo,hi,_)in REGIONS.items():
                        x.setRange(reg,lo,hi)
                        actual=full.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.NormSet(ROOT.RooArgSet(x)),ROOT.RooFit.Range(reg)).getVal()
                        analytic=w.function(f'local_acceptance_{proc}_{era}_{reg}').getVal()
                        assert abs(analytic-actual)<1e-8
                        fractions.append(analytic)
                        mu_checks.append(dict(era=era,process=proc,region=reg,mean=mv.getVal(),difference=analytic-actual,change=analytic-original_fractions[reg]))
                    assert abs(sum(fractions)-1)<1e-10
                    assert any(abs(fractions[k]-original_fractions[reg])>1e-7 for k,reg in enumerate(REGIONS))
                mv.setVal(central)
    if deep:
        def digest(data,remap):
            rows=[]
            for k in range(data.numEntries()):
                row=data.get(k);ch=row.find('CMS_channel').getLabel()
                if remap:
                    reg,channel=ch.split('_',1);xx=row.getRealValue('mass_'+reg)
                else:channel=ch;xx=row.getRealValue('mass')
                rows.append((channel,xx,data.weight()))
            h=hashlib.sha256()
            for channel,xx,weight in sorted(rows):h.update(channel.encode()+struct.pack('dd',xx,weight))
            return h.hexdigest()
        assert digest(w.data('data_obs'),True)==digest(old.data('data_obs'),False)
    source.Close();handle.Close();assert sha(path)==before
    return dict(mass=float(mass),passed=True,workspace_sha256=before,means_floating=8,
                original_mean_bounds_preserved=True,original_widths_tails_fixed=True,
                discrete_indices=4,twice_NLL_penalty=32.,bins=1800,events=318924,
                partition_max_error=max(rate_partitions),deep_audit=deep,
                intensity_max_relative_error=max(intensity_errors,default=0),dynamic_mu_checks=mu_checks,
                source_sha256=meta['source_sha256'],source_unchanged=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--masses',nargs='+');parser.add_argument('--all',action='store_true');args=parser.parse_args()
    masses=[f'{k/10:.1f}'for k in range(60,121)]if args.all else args.masses or ['7.5','9.0']
    rows=[audit(m,deep=m in ['7.5','9.0'])for m in masses]
    result=dict(passed=all(q['passed']for q in rows),workspaces=len(rows),new_fits=0,rows=rows,
                interpretation='Prior CR estimates initialize floating means only; no CR measurement constraint. Degree7 global intensity and two-region normalization audited.',
                builder_sha256=sha(Path(__file__).with_name('model.py')))
    path=R/('all_mass_model_audit.json'if args.all else 'pilot_model_delta_audit.json')
    path.write_text(json.dumps(result,indent=2)+'\n');print(path,flush=True)

"""Attach nominal-fit diagnostics using a parity-validated private Combine build."""
from pathlib import Path
import argparse,json,time,hashlib,traceback
import ROOT
from run_observed import R,P,execute,extract,atomic
from gof_diagnostics import read_nominal_diagnostics,DIAGNOSTIC_COMBINE,RTD_FLAG
from observed_helpers import dump
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')

def audit(path):
    D=path.parent;J=D/'joint';out=D/'gof_diagnostics.json'
    if out.exists():return
    data=json.loads(path.read_text())
    original=J/'gof_observed';receipt=original/'combine.receipt.json'
    if not receipt.exists():return
    original_receipt=json.loads(receipt.read_text())
    if original_receipt['exit_code']!=0 or ('gof_observed' not in data and (D/'completion.json').exists()):
        dump(out,dict(mass=data['mass'],audit_completed=True,available=False,valid=False,
             flags=['original_GoF_failed_or_missing_score'],original_receipt=original_receipt))
        extract(f'{data["mass"]:.1f}');return
    if 'gof_observed' not in data:return
    cmd=json.loads((original/'combine.command.json').read_text())
    source=Path(cmd[cmd.index('-d')+1]);before=hashlib.sha256(source.read_bytes()).hexdigest()
    cmd[0]=str(DIAGNOSTIC_COMBINE);cmd+=['--X-rtd',RTD_FLAG]
    directory=J/'gof_observed_diagnostic';done=execute(cmd,directory,timeout=1200)
    files=sorted(directory.glob('higgsCombine*.root'))
    result=read_nominal_diagnostics(files[-1],source,data['gof_observed']) if files else dict(available=False,valid=False,flags=['diagnostic_output_missing'])
    result.update(mass=data['mass'],audit_completed=True,receipt=done,source=str(source),source_sha256=before,
                  source_unchanged=hashlib.sha256(source.read_bytes()).hexdigest()==before)
    fd=data.get('fit',{}).get('s',{})
    if fd.get('point_valid') and result.get('deviance_with_aux') is not None:
        result['FD_reference_deviance']=fd['exact_deviance']
        result['nominal_minus_FD_deviance']=result['deviance_with_aux']-fd['exact_deviance']
        bound=float(cmd[cmd.index('--rMax')+1])
        result['FD_reference_within_GoF_range']=abs(fd['r'])<bound
        if result['FD_reference_within_GoF_range'] and result['nominal_minus_FD_deviance']>.01:
            result['flags'].append('known_better_FitDiagnostics_minimum')
            result['valid']=False
    dump(out,result);extract(f'{data["mass"]:.1f}')
    print('GoF nominal audit',data['mass'],result.get('valid'),result.get('flags'),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');a=p.parse_args()
    while True:
        files=sorted((R/'observed').glob('m*/results.json'))
        for path in files:
            try:audit(path)
            except Exception as error:
                failure=dict(error=str(error),traceback=traceback.format_exc(),audit_completed=True,
                             available=False,valid=False,flags=['diagnostic_audit_exception'])
                atomic(path.parent/'gof_diagnostic_failure.json',failure)
                dump(path.parent/'gof_diagnostics.json',failure)
        if not a.watch:break
        if len(list((R/'observed').glob('m*/gof_diagnostics.json')))==61:break
        time.sleep(30)

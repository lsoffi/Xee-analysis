"""Read-only native-limit minimum audit for shared SR7/CR4 with floating Upsilon means."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import re
import ROOT

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'
C=Path('/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/CombinedLimit')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    rows=[]
    for receipt in sorted(R.glob('observed/m*/joint/limits/combine.receipt.json')):
        directory=receipt.parent
        q=json.loads((directory.parent.parent/'results.json').read_text())
        command=json.loads((directory/'combine.command.json').read_text())
        log=(directory/'combine.log').read_text(errors='replace')
        native=json.loads(receipt.read_text())
        data=re.search(r'NLL at global minimum of data = ([^ ]+) \(r=([^\)]+)\)',log)
        asimov=re.search(r'NLL at global minimum of asimov = ([^ ]+) \(r=([^\)]+)\)',log)
        files=list(directory.glob('higgsCombine*.root'))
        schema={}
        if files:
            f=ROOT.TFile.Open(str(files[-1]));tree=f.Get('limit')
            schema=dict(keys=[key.GetName() for key in f.GetListOfKeys()],branches=[b.GetName() for b in tree.GetListOfBranches()],entries=int(tree.GetEntries()))
            f.Close()
        s=q['fit']['s'];b=q['fit']['b'];rhat=float(data[2]) if data else None
        qt='--qtilde' not in command or command[command.index('--qtilde')+1].lower() not in ['0','false']
        if not s['point_valid']:
            reference_kind='FD_signed_minimum_not_numerically_valid'
            reference_r=None;reference_D=None
        elif s['r']>=0:
            reference_kind='valid_positive_FD_signal_plus_background'
            reference_r=s['r'];reference_D=s['exact_deviance']
        elif qt and b['point_valid']:
            reference_kind='physical_boundary_zero_uses_valid_FD_background_reference'
            reference_r=0.;reference_D=b['exact_deviance']
        else:
            reference_kind='no_usable_physical_reference'
            reference_r=None;reference_D=None
        flags=[]
        if data is None:flags.append('global_data_minimum_not_logged')
        if reference_r is None:flags.append('no_valid_FD_reference_for_minimum_validation')
        delta=rhat-reference_r if rhat is not None and reference_r is not None else None
        normalized=abs(delta)/s['error'] if delta is not None and s['r']>0 and s.get('error') and s['error']>0 else None
        rows.append(dict(mass=q['mass'],exit_code=native['exit_code'],qtilde=qt,
            Asymptotic_rhat=rhat,Asymptotic_logged_NLL=data[1] if data else None,
            Asimov_rhat=float(asimov[2]) if asimov else None,Asimov_logged_NLL=asimov[1] if asimov else None,
            FD_signed_rhat=s['r'],FD_point_valid=s['point_valid'],FD_EDM=s['edm'],FD_error=s.get('error'),
            reference_kind=reference_kind,physical_reference_r=reference_r,physical_reference_deviance=reference_D,
            delta_r=delta,abs_delta_r_over_FD_Hesse_error=normalized,
            observed_CLs_crossing_verified=q.get('observed_limit_check',{}).get('verified',False),
            unconditional_minimum_fully_verified=False,unconditional_parameter_snapshot_available=False,
            unconditional_EDM_status_available=False,absolute_NLL_comparison_performed=False,
            flags=flags,root_schema=schema,
            files={str(path):sha(path) for path in [receipt,directory/'combine.command.json',directory/'combine.log',*files]}))
    rows.sort(key=lambda q:q['mass'])
    positive=[q for q in rows if q['reference_kind']=='valid_positive_FD_signal_plus_background']
    boundary=[q for q in rows if q['reference_kind']=='physical_boundary_zero_uses_valid_FD_background_reference']
    output=dict(audit_date=datetime.datetime.now().astimezone().isoformat(),read_only=True,new_fits=0,
                model_configuration='Shared SR Chebyshev7/Bernstein7 envelope + CR Chebyshev4; eight Upsilon means floating',
                new_likelihood_minimizations=0,new_limit_estimates=0,production_modified=False,
                completed_limit_files=len(rows),positive_FD_references=len(positive),physical_boundary_references=len(boundary),
                no_valid_FD_reference=[q['mass'] for q in rows if q['physical_reference_r'] is None],
                maximum_positive_r_difference_over_Hesse_error=max((q['abs_delta_r_over_FD_Hesse_error'] for q in positive if q['abs_delta_r_over_FD_Hesse_error'] is not None),default=None),
                maximum_boundary_rhat=max((q['Asymptotic_rhat'] for q in boundary if q['Asymptotic_rhat'] is not None),default=None),
                conclusion=f'Logged rhat can be compared with a usable physical FD reference at {len(positive)+len(boundary)} of {len(rows)} completed points. The differences are tabulated without a new quality gate. This partial check does not prove the same global minimum: full saved parameters, indices, EDM and unshifted NLL are absent.',
                source_findings=[
                    dict(file=str(C/'src/AsymptoticLimits.cc'),line=166,finding='The native observed fit resets r to 0.1*rMax and restricts r>=0 by default qtilde.'),
                    dict(file=str(C/'src/AsymptoticLimits.cc'),line=190,finding='The observed global minimum is fitted, then parameters are stored only in an in-memory CheapValueSnapshot. The minimize return is not retained or checked here.'),
                    dict(file=str(C/'src/AsymptoticLimits.cc'),line=199,finding='Only NLL and r are logged with %g, i.e. default six significant digits.'),
                    dict(file=str(C/'src/AsymptoticLimits.cc'),line=205,finding='Full fitFreeD parameter printing requires verbose>1; current commands use verbose1.'),
                    dict(file=str(C/'interface/AsymptoticLimits.h'),line=65,finding='fitFreeD/fitFreeA are CheapValueSnapshot members; no persisted RooFitResult for the global fit.'),
                    dict(file=str(C/'src/CachingNLL.cc'),line=383,finding='CachingAddNLL initializes an arbitrary constantZeroPoint from its initial state; the printed NLL is not directly comparable across methods.'),
                    dict(file=str(C/'src/CachingNLL.cc'),line=653,finding='That constant offset contributes to evaluate() unless explicitly disabled.'),
                    dict(file=str(C/'src/AsymptoticLimits.cc'),line=311,finding='CLs uses differences from the native minNllD. A correct CLs crossing alone cannot verify that denominator minimum.'),
                ],
                limitations=[
                    'rhat agreement does not establish agreement of all nuisance parameters, family indices, NLL or numerical convergence.',
                    'A signed FD estimate below zero is not a discrepancy with a qtilde limit at r=0; the valid background fit is the appropriate known physical reference.',
                    'The rounded, offset NLL cannot be subtracted from FD.minNll or from the exact absolute deviance without reconstructing the identical native zero point.',
                    'Generic minimizer status lines in combine_logger.out are interleaved across expected-quantile scans and conditional fits. The current verbosity does not identify a complete saved convergence record for the observed global fit; the last preceding generic status must not be treated as that record.',
                    'No claim is made that the current limits are wrong; minimum validation remains incomplete. Native limits and expected bands were not altered.',
                ],rows=rows,
                source_sha256={str(C/name):sha(C/name) for name in ['src/AsymptoticLimits.cc','interface/AsymptoticLimits.h','src/CachingNLL.cc']})
    (R/'limits_unconditional_minimum_audit.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    max_positive=output['maximum_positive_r_difference_over_Hesse_error']
    max_boundary=output['maximum_boundary_rhat']
    # Early partial scans may contain only one physical-reference class.
    fmt=lambda value: f'{value:.6g}' if value is not None else 'non disponibile'
    unavailable=', '.join(f'{mass:.1f}' for mass in output['no_valid_FD_reference']) or 'nessuna'
    paragraphs=[
        '# Audit del minimo incondizionato: SR7/CR4, medie Υ libere',
        '',
        f'Aggiornamento: {output["audit_date"]}. Controllati {len(rows)} output AsymptoticLimits già prodotti. Nessun fit, rivalutazione dei limiti, modifica delle bande o nuovo criterio di selezione.',
        '',
        '## Cosa è verificabile',
        '',
        f'I log salvano r̂ del fit osservato incondizionato. Per {len(positive)} punti il riferimento FD è positivo e numericamente valido; lo scarto massimo è {fmt(max_positive)} volte il suo errore Hesse. È un confronto numerico descrittivo, non un nuovo test statistico.',
        '',
        f'Per altri {len(boundary)} punti il FD firmato è negativo. AsymptoticLimits usa q̃ con r≥0: il riferimento ammissibile è quindi il fit B-only a r=0. Il massimo r̂ nativo in questi casi è {fmt(max_boundary)}. Non sono stati aggiunti flag per la differenza rispetto al FD negativo.',
        '',
        f'Masse senza un riferimento FD s+b numericamente valido: {unavailable} GeV. Un r̂ diverso in questi punti non dimostra, da solo, che AsymptoticLimits abbia perso un minimo valido.',
        '',
        '## Cosa rimane non verificato',
        '',
        'L’accordo in r̂ non dimostra che siano uguali le nuisance, la famiglia selezionata o la NLL. Gli output ROOT contengono la tabella dei limiti e dei quantili, senza RooFitResult, snapshot del minimo incondizionato, indici di famiglia, stato o EDM. Il codice conserva i parametri in un CheapValueSnapshot soltanto in memoria; la loro stampa richiede verbosità maggiore di quella usata.',
        '',
        'La NLL stampata usa sei cifre significative e include l’offset costante inizializzato da CachingAddNLL. Non è quindi confrontabile direttamente con FD.minNll o con la devianza assoluta. Anche le righe generiche di stato del minimizzatore appartengono a numerosi fit interni e non forniscono un record completo identificato del minimo osservato.',
        '',
        'Il crossing CLs≈0,05 verifica la ricerca del limite rispetto al denominatore nativo usato; non prova che quel denominatore corrisponda al minimo globale. Gli scarti di r̂ sono riportati descrittivamente; la verifica completa NLL/stato rimane indisponibile dai soli output conservati. Questo audit non conclude che i limiti siano errati.',
        '',
        'Il JSON conserva i valori per massa, lo schema ROOT, gli hash dei file e i riferimenti alle righe del codice. Il parser `scripts/shared_sr7_floatmu_fitlimits/audit_limits_minimum_readonly.py` può essere rieseguito sui 61 risultati finali: legge soltanto i prodotti esistenti e riscrive questi due documenti di audit.',
        '',
        '[Dettaglio numerico e riferimenti](limits_unconditional_minimum_audit.json)',
    ]
    (R/'limits_unconditional_minimum_audit.md').write_text('\n'.join(paragraphs)+'\n')
    print(json.dumps({k:output[k] for k in ['completed_limit_files','positive_FD_references','physical_boundary_references','no_valid_FD_reference','maximum_positive_r_difference_over_Hesse_error','maximum_boundary_rhat']},indent=2))


if __name__=='__main__':main()

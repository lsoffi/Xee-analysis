import pathlib,json,subprocess,time,csv,hashlib
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_background_validation';S=P/'scripts/validation';plotpython='/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/fit-env/bin/python'
while True:
 try:
  a=json.loads((R/'postprocess_status.json').read_text());b=json.loads((R/'retry_progress.json').read_text());ready=a['completed_fit_batches']==200 and b['baseline_complete'] and b['active']==0 and b['pending_unsubmitted']==0
 except (FileNotFoundError,json.JSONDecodeError):ready=False
 if ready:break
 time.sleep(10)
for interpreter,script in [(str(P/'bin/python'),'audit_toys.py'),(plotpython,'summarize_toys.py'),(plotpython,'summarize_effective.py'),(plotpython,'plot_validation_summary.py'),(plotpython,'build_report.py')]:
 subprocess.run([interpreter,str(S/script)],check=True)
rows=list(csv.DictReader((R/'toys_all.csv').open()));effective=list(csv.DictReader((R/'toys_effective.csv').open()));assert len(rows)==len(effective)==2000
assert len({(r['mass'],r['generator'],r['injection'],r['seed'],r['toy']) for r in rows})==1000
receipts=list(R.glob('toys/m*/*/*/batch_*/*/receipt.json'));assert len(receipts)==300
assert all(json.loads(p.read_text())['exit_code']==0 for p in receipts)
audit=json.loads((R/'model_integrity_audit.json').read_text());verified=[]
for m in ['7.5','9.5']:
 src=P/f'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/{m}/Xee_ee_4_allYears.root';sha=hashlib.sha256(src.read_bytes()).hexdigest();assert all(r['input_sha256']==sha for r in audit if r['mass']==m);verified.append(dict(mass=m,source=str(src),sha256=sha,unchanged=True))
retries=[json.loads(p.read_text()) for p in (R/'retries').rglob('receipt.json')];status=dict(datasets=1000,fit_results=2000,initial_fit_batches=200,generation_batches=100,initial_execution_failures=0,initial_usable_intervals=sum(r.get('valid')=='True' for r in rows),effective_usable_intervals=sum(r.get('valid')=='True' for r in effective),effective_usable_point_fits=sum(r.get('valid_point')=='True' for r in effective),retry_attempts=len(retries),retry_timeouts=sum(r['exit_code']==124 for r in retries),retry_interrupted=sum(r['exit_code']==-2 for r in retries),input_files=verified,report=str(R/'REPORT.md'),completed=True)
(R/'completion.json').write_text(json.dumps(status,indent=2));(P/'runs/LATEST_BACKGROUND_VALIDATION.txt').write_text(str(R)+'\n');print(json.dumps(status,indent=2),flush=True)

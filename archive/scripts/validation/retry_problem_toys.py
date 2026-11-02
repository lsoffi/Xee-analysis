import pathlib,csv,json,subprocess,time,concurrent.futures
R=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation');start=time.time();seen=set();active={};completed=0

def retry(row):
 mass=str(row['mass']);g=row['generator'];inj='zero' if float(row['injection'])==0 else ('moderate' if float(row['injection'])==3 else 'signal');seed=int(row['seed']);toy=int(row['toy']);model=row['model']
 src=next((R/f'toys/m{mass}/{g}/{inj}').glob(f'batch_*/{model}/command.json')) if False else None
 for path in (R/f'toys/m{mass}/{g}/{inj}').glob('batch_*/scenario.json'):
  if json.loads(path.read_text())['seed']==seed:src=path.parent/model/'command.json';break
 assert src
 w=R/f'retries/m{mass}/{g}/{inj}/{model}/seed{seed}_toy{toy:02d}';w.mkdir(parents=True,exist_ok=True)
 if (w/'receipt.json').exists():
  old=json.loads((w/'receipt.json').read_text())
  if old['exit_code']!=-2:return str(w)
  (w/'interrupted_receipt.json').write_text(json.dumps(old,indent=2))
 cmd=json.loads(src.read_text());cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='1';cmd[cmd.index('-n')+1]='_retry';cmd+=['--robustFit','1','--cminDefaultMinimizerTolerance','0.001'];tf=cmd[cmd.index('--toysFile')+1];single=w/'single_toy.root';subprocess.run(['/Users/imb/Documents/Codex/Analisi-statistica-Xee/bin/python','/Users/imb/Documents/Codex/Analisi-statistica-Xee/scripts/validation/extract_single_toy.py',tf,str(toy),str(single)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);cmd[cmd.index('--toysFile')+1]=str(single);cmd[cmd.index('-t')+1]='1'
 (w/'command.json').write_text(json.dumps(cmd,indent=2));(w/'initial_result.json').write_text(json.dumps(row,indent=2));t=time.time()
 with (w/'combine.log').open('w') as f:
  try:rc=subprocess.run(cmd,cwd=w,stdout=f,stderr=subprocess.STDOUT,timeout=180).returncode
  except subprocess.TimeoutExpired:rc=124
 (w/'receipt.json').write_text(json.dumps(dict(exit_code=rc,elapsed_seconds=time.time()-t),indent=2));return str(w)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
 while True:
  for fu,key in list(active.items()):
   if fu.done():
    try:fu.result()
    except Exception as e:print(str(e),flush=True)
    completed+=1;del active[fu]
  status=json.loads((R/'postprocess_status.json').read_text()) if (R/'postprocess_status.json').exists() else {};full=status.get('completed_fit_batches',0)>=200;limit=8 if full else 2
  rows=list(csv.DictReader((R/'toys_all.csv').open())) if (R/'toys_all.csv').exists() else []
  pending=[]
  for r in rows:
   if r.get('valid')!='False' or 'toy' not in r:continue
   if float(r['injection'])==25 and int(r['seed']) not in [300510,300530]:continue
   if r['model']=='original' and (int(r['b_covQual'])!=3 or int(r['s_covQual'])!=3):continue
   if int(r['seed'])%10!=0:continue
   key=(r['mass'],r['generator'],r['injection'],r['seed'],r['toy'],r['model'])
   if key not in seen:pending.append((key,r))
  for key,row in pending[:max(0,limit-len(active))]:seen.add(key);active[pool.submit(retry,row)]=key
  (R/'retry_progress.json').write_text(json.dumps(dict(completed=completed,active=len(active),pending_unsubmitted=max(0,len(pending)-max(0,limit-len(active))),baseline_complete=full,elapsed=time.time()-start),indent=2))
  if full and not active and not pending:break
  time.sleep(3)

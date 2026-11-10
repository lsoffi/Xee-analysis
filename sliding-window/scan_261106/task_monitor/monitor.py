import pathlib,subprocess,json,time,datetime,html,os
base=pathlib.Path('/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106')
out=pathlib.Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/task_monitor');out.mkdir(parents=True,exist_ok=True)
local=base/'task_monitor'
import fcntl
lock=open(local/'monitor.lock','w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit('Monitor already active')
(local/'pid').write_text(str(os.getpid()))
while True:
 stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();data={'updated_utc':stamp,'tasks':[]}
 for name,query,pattern,target,title in [('gof_sr_single_toy','ClusterId==13692094','results/*.json',8000,'GoF SR phi–J/psi'),('high_c4l4_diagnostics','ClusterId==13692091','m*/diagnostics.json',50,'FitDiagnostics C4+L4 4.1–9.0 GeV')]:
  task={'name':title,'directory':name,'expected':target,'warnings':[]};p=base/name
  try:
   r=subprocess.run(['condor_q','-name','bigbird18.cern.ch','-constraint',query,'-json'],capture_output=True,text=True,timeout=45)
   if r.returncode:raise RuntimeError(r.stderr[:300])
   jobs=json.loads(r.stdout or '[]');counts={}
   for j in jobs:
    status=str(j.get('JobStatus'));counts[status]=counts.get(status,0)+1
    if j.get('JobStatus')==5:task['warnings'].append('HOLD: '+j.get('HoldReason','unknown'))
   task['queue_status']=counts;task['cpu_seconds']=sum(j.get('RemoteUserCpu',0)+j.get('RemoteSysCpu',0) for j in jobs)
  except Exception as e:task['warnings'].append('Queue unavailable: '+str(e));jobs=None
  results=list(p.glob(pattern));task['completed']=len(results)
  if name=='gof_sr_single_toy':
   groups={}
   for f in results:
    try:
     d=json.load(open(f));groups.setdefault(d['mass'],[]).append(d['statistic'])
    except Exception as e:task['warnings'].append(str(e))
   task['pvalues']=[]
   for m,ts in sorted(groups.items()):
    obs=json.load(open(base/'gof_sr_only'/'results'/f'{m:.1f}_observed.json'))['statistics'][0]
    task['pvalues'].append({'mass':m,'toys':len(ts),'tail':sum(t>=obs for t in ts),'p':sum(t>=obs for t in ts)/len(ts),'preliminary':len(ts)<500})
   exits=list((p/'logs').glob('*.exit'))
   failed=[f.name for f in exits if f.read_text().strip()!='0']
   if failed:task['warnings'].append(str(len(failed))+' nonzero exit receipts; examples: '+str(failed[:10]))

  else:
   task['workspaces']=len(list(p.glob('m*/envelope.root')));task['fit_files']=len(list(p.glob('m*/fitDiagnostics*.root')))
   if jobs is not None:
    accounted=len(results)+len(jobs)
    if accounted<target:task['warnings'].append(str(target-accounted)+' points without active job or completed report; inspect Condor history for failures.')
   for f in results:
    try:
     d=json.load(open(f))
     for s,v in d['stages'].items():
      if v['status']!=0:task['warnings'].append(f'm{d["mass"]} fit_{s} status {v["status"]}')
    except Exception as e:task['warnings'].append(str(e))
  if jobs==[] and len(results)<target:task['warnings'].append('No active jobs but results incomplete. Requires diagnosis; no automatic resubmission.')
  for f in p.glob('logs/*.err'):
   try:
    s=f.read_text(errors='replace')
    if 'Traceback' in s or 'Segmentation violation' in s:task['warnings'].append('Error in '+f.name)
   except OSError:pass
  data['tasks'].append(task)
 p=base/'high_c4l4_limits_bias'
 extra={'name':'High-mass limits and C4/L4 bias','directory':'high_c4l4_limits_bias','expected':40000,'completed':len(list(p.glob('m*/results/*.json'))),'limits':len(list((p/'limits').glob('*.json'))),'generator_snapshots':len(list(p.glob('m*/gen_*.root'))),'warnings':[]}
 try:
  q=subprocess.run(['condor_q','-name','bigbird18.cern.ch','-constraint','DAGManJobId==13692098','-json'],capture_output=True,text=True,timeout=45)
  if q.returncode:raise RuntimeError(q.stderr[:200])
  js=json.loads(q.stdout or '[]');extra['queue_status']={}
  for j in js:
   k=str(j.get('JobStatus'));extra['queue_status'][k]=extra['queue_status'].get(k,0)+1
   if j.get('JobStatus')==5:extra['warnings'].append(j.get('HoldReason','held'))
 except Exception as e:extra['warnings'].append(str(e))
 data['tasks'].append(extra)
 p=base/'upsilon_limit_diagnostics'
 extra={'name':'Upsilon observed-limit diagnosis 8.6 / 9.0','expected':108,'completed':len(list(p.glob('results/*/result.json'))),'warnings':[]}
 try:
  q=subprocess.run(['condor_q','-name','bigbird18.cern.ch','-constraint','ClusterId==13692444','-json'],capture_output=True,text=True,timeout=45)
  if q.returncode:raise RuntimeError(q.stderr[:200])
  extra['queue_status']={}
  for j in json.loads(q.stdout or '[]'):
   k=str(j.get('JobStatus'));extra['queue_status'][k]=extra['queue_status'].get(k,0)+1
   if j.get('JobStatus')==5:extra['warnings'].append(j.get('HoldReason','held'))
  bad=[str(f.parent.name) for f in p.glob('results/*/exit_code.txt') if f.read_text().strip()!='0']
  if bad:extra['warnings'].append('Failed jobs: '+str(bad))
 except Exception as e:extra['warnings'].append(str(e))
 data['tasks'].append(extra)
 text=json.dumps(data,indent=2)
 for dest in [local/'status.json',out/'status.json']:
  tmp=dest.with_suffix('.tmp');tmp.write_text(text);tmp.replace(dest)
 page='<html><head><meta http-equiv="refresh" content="300"></head><body><h1>Xee task monitor</h1><p>Updated UTC: '+stamp+'</p><p>Checks every 5 minutes. Queue states: 1 idle, 2 running, 3 removed, 4 completed, 5 held. This monitor reports issues; it does not restart jobs.</p>'
 for t in data['tasks']:page+='<h2>'+html.escape(t['name'])+'</h2><pre>'+html.escape(json.dumps(t,indent=2))+'</pre>'
 page+='</body></html>';tmp=out/'index.tmp';tmp.write_text(page);tmp.replace(out/'index.html')
 print(stamp,[(t['directory'],t['completed'],t.get('queue_status'),t['warnings']) for t in data['tasks']],flush=True)
 if all(t['completed']==t['expected'] and not t.get('queue_status') and not t['warnings'] for t in data['tasks']):break
 time.sleep(300)

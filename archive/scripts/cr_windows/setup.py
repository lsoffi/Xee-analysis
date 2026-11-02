from pathlib import Path
import math,json
P=Path.cwd();S=P/'scripts/cr_windows';R=P/'runs/20261029_cr_windows';configs=[]
for tag,lo,hi in [('cr90_106',9.,10.6),('cr91_105',9.1,10.5),('cr92_104',9.2,10.4)]:
 il=math.floor((lo-6)*75+1e-8);ih=math.ceil((hi-6)*75-1e-8);low=6+il/75;high=6+ih/75;nb=ih-il;configs.append(dict(tag=tag,requested=[lo,hi],range=[low,high],bins=nb))
 for degree in [2,3]:
  src=P/'scripts'/('upsilon_control' if degree==2 else 'upsilon_control_cheb3');ss=S/tag/f'degree{degree}';ss.mkdir(parents=True,exist_ok=True);rr=R/tag/f'degree{degree}';rr.mkdir(parents=True,exist_ok=True)
  for name in ['prepare_control.py','build_joint.py','run_joint.py','extract.py','report.py']:
   if degree==2 and name!='prepare_control.py':continue
   s=(src/name).read_text().replace('20261029_upsilon_control_cheb3',f'20261029_cr_windows/{tag}/degree{degree}') if degree==3 else (src/name).read_text().replace('20261029_upsilon_control',f'20261029_cr_windows/{tag}/degree{degree}')
   if name=='prepare_control.py':s=s.replace('win=dict(low95_common_GeV=8.8,high95_common_GeV=10.8);LO=8.8;HI=10.8;NB=150',f'win=dict(low95_common_GeV={lo},high95_common_GeV={hi});LO={low!r};HI={high!r};NB={nb}').replace('CONTROL 8.8-10.8',f'CONTROL {low}-{high}')
   if name=='build_joint.py':s=s.replace("R/'signal/workspace.root'","P/'runs/20261029_upsilon_control_cheb3/signal/workspace.root'").replace('9.8,8.8,10.8',f'9.8,{low!r},{high!r}').replace('crmass.setBins(150)',f'crmass.setBins({nb})').replace('CR_range=[8.8,10.8]',f'CR_range=[{low!r},{high!r}]')
   if name=='report.py':
    s=s.replace("P/'runs/20261029_upsilon_control/results.json'","P/'runs/20261029_upsilon_control_cheb3/results.json'").replace('Grado 2 + controllo Υ','CR 8,8–10,8').replace('Grado 3 + controllo Υ',f'CR {low:.3f}–{high:.3f}').replace('confronto grado 2 / grado 3','confronto finestre CR').replace('(8.8-r[',f'({low!r}-r[').replace('(10.8-r[',f'({high!r}-r[').replace('CR8,8–10,8',f'CR{low:.4f}–{high:.4f}')
   (ss/name).write_text(s)
(R/'configs.json').write_text(json.dumps(configs,indent=2))

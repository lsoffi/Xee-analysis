from pathlib import Path
import subprocess
old=Path('low_logpoly3');new=Path('low_logpoly3_totalC0');new.mkdir();(new/'logs').mkdir()
s=(old/'build.py').read_text()
key="            amplitude=form('logpoly_anchor_'+era,f'{density:.17g}*@0',[d['node']])"
assert s.count(key)==1
s=s.replace(key,key+'''
            tails=[]
            for res in resonances:
                if res['group']==g:continue
                t=f'(({edge:.17g}-@0)/abs(@1))'
                left=f'exp(-0.5*@2*@2-@3*log(1+@2*(-@2-{t})/@3))'
                right=f'exp(-0.5*@4*@4-@5*log(1+@4*({t}-@4)/@5))'
                core=f'exp(-0.5*{t}*{t})'
                pdf=f'(({t}<-@2)?{left}:(({t}>@4)?{right}:{core}))'
                tails.append(form('edge_tail_'+res['name']+'_'+era,
                    '@6/(abs(@1)*@7)*'+pdf,res['params']+[res['norm'],res['den']]))
            excess=add('edge_extra_density_'+era,tails)
            d['node']=form('N_cont_totalC0_'+g+'_'+era,f'@0+@1/({density:.17g})',[d['node'],excess])
            for entry in source_meta:
                if entry['group']==g:entry['continuum_norm']=d['node'].GetName()
''')
s=s.replace("bg=var('N_cont_'+group+'_'+era","bg=var(('anchor_scale_' if group=='omega_phi' else 'N_cont_')+group+'_'+era")
s=s.replace('metadata=dict(degree','metadata=dict(continuity="total physical event density",degree')
(new/'build.py').write_text(s)
for name in ['worker.sh','draw.py','jobs.sub']:
 s=(old/name).read_text().replace('low_logpoly3',new.name).replace('C0 ai bordi comuni','C0 totale al bordo')
 (new/name).write_text(s)
(new/'worker.sh').chmod(0o755)
(new/'masses.txt').write_text('0.6\n')
subprocess.run(['condor_submit','jobs.sub'],cwd=new,check=True)

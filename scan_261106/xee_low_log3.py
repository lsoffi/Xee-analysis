from pathlib import Path
import json,subprocess
base=Path.cwd()
p=base/'low_logpoly3';p.mkdir()
(p/'logs').mkdir()
rec=[r for r in json.loads((base/'records.json').read_text()) if r['mass']<.65]
assert [r['mass'] for r in rec]==[.4,.5,.6]
(p/'masses.txt').write_text('0.4\n0.5\n0.6\n')
s=(base/'build.py').read_text().replace('BASE=Path(__file__).resolve().parent','BASE=Path(__file__).resolve().parent.parent')
a=s.index('        matrix=[]');b=s.index('        # Original signal',a)
s=s[:a]+'''        assert len(contacts)<=1
        rank=len(contacts); targets=[]; C=np.empty((0,4))
        anchor=0.
        if contacts:
            g,edge=contacts[0];d=norms[g]
            anchor=2*(edge-a)/width-1
            density=math.exp(d['slope']*(edge-d['lo']))/d['integral']
            targets=[(g,density)]
            amplitude=form('logpoly_anchor_'+era,f'{density:.17g}*@0',[d['node']])
        else:
            scale=sum(float(z['events']) for z in srrows)/width
            amplitude=var('logpoly_density_'+era,scale,scale*1e-6,scale*100)
        shape=[var(f'logpoly_{era}_c{j}',0,-8,8,.02) for j in (1,2,3)]
        coeff=[amplitude]+shape
        free=shape if contacts else coeff
        nodes,weights=np.polynomial.legendre.leggauss(8)
'''+s[b:]
a=s.index('                    weights=[]',s.index('        for reg,rows,group'));b=s.index("                processes.append(('continuum'",a)
s=s[:a]+'''                    terms=[]
                    for u,wgt in zip(nodes,weights):
                        t=2*((lo+hi)/2+(hi-lo)*u/2-a)/width-1
                        exponent='+'.join(f'({t**j-anchor**j:.17g})*@{j}' for j in (1,2,3))
                        terms.append(f'({wgt*(hi-lo)/2:.17g})*exp({exponent})')
                    parts.append(form('contbin_'+label+'_'+str(k),'@0*('+'+'.join(terms)+')',coeff))
'''+s[b:]
s=s.replace('exact_bin_integration=True','SR_model="exp(cubic in scaled mass)",SR_quadrature_order=8,exact_bin_integration=False')
(p/'build.py').write_text(s)
draw=(base/'early_fd/draw_one.py').read_text().replace("scan_SR_CR/fit'","scan_SR_CR/low_logpoly3'").replace('SR Chebyshev 6','SR exp(pol. grado 3)')
(p/'draw.py').write_text(draw)
s=(base/'worker.sh').read_text().replace('XEE_DEGREE=6','XEE_DEGREE=3').replace('out="$run/results','out="$run/low_logpoly3/results').replace('"$run/build.py"','"$run/low_logpoly3/build.py"').replace('in Significance AsymptoticLimits FitDiagnostics','in FitDiagnostics').replace('exit "$rc"','''cp model.json workspace.root "$out/"
if [[ $rc == 0 ]]; then
 cd "$run"
 python3 low_logpoly3/draw.py "$m" "$out"
fi
exit "$rc"''')
(p/'worker.sh').write_text(s);(p/'worker.sh').chmod(0o755)
(p/'jobs.sub').write_text((base/'jobs.sub').read_text())
subprocess.run(['condor_submit','jobs.sub'],cwd=p,check=True)

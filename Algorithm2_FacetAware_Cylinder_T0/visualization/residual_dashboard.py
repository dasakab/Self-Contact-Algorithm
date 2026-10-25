#!/usr/bin/env python3
"""Live cylinder FEM residual dashboard. Standard library only.

 python residual_dashboard.py --follow retry31.log
 python residual_dashboard.py --replay retry31.log

Creates monitor.html (browser refreshes every 8 seconds), and SVG/CSV
for each Newton iteration curve. Parses only actual solver log lines.
Does not run the Fortran solver or modify its calculation.
"""
import argparse, collections, csv, html, math, os, pathlib, re, sys, time
from datetime import datetime

INC = re.compile(r'INC NO\s*=\s*(\d+)\s+ITR NO\s*=\s*(\d+)')
RES = re.compile(r'UNBALANCE FORCE\s*=\s*([-+0-9.eEdD]+)')
TOL = re.compile(r'TOLEC\s*=\s*([-+0-9.eEdD]+)')
SC  = re.compile(r'SELF_CONTACT_STATUS inc inner outer:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
RING= re.compile(r'SC_RING_(IN|OUT)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+0-9.eEdD]+)')

def f(s):
    return float(s.replace('D','E').replace('d','e'))

class Watch:
    def __init__(self, workdir):
        self.workdir=pathlib.Path(workdir); self.cur=None; self.iter=None; self.tol=5e-4
        self.data=collections.defaultdict(dict); self.conv={}; self.contact={}; self.rings=collections.defaultdict(set)
        self.stopped=False;self.stopped_reason='';self.last_line='';self.changed=False;self.events=[];self.last_render=0
    def event(self,message):
        stamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        entry=f'[{stamp}] {message}'
        self.events.append(entry);self.events=self.events[-60:]
        print(entry,flush=True); self.changed=True
    def process(self,line):
        m=INC.search(line)
        if m:
            self.cur,self.iter=map(int,m.groups());self.changed=True;return
        m=RES.search(line)
        if m and self.cur is not None and self.iter is not None:
            self.data[self.cur][self.iter]=f(m.group(1)); self.changed=True;return
        m=TOL.search(line)
        if m:self.tol=f(m.group(1));return
        m=SC.search(line)
        if m:
            inc,ia,ih,oa,oh=map(int,m.groups());was=self.contact.get(inc,(0,0,0,0))
            self.contact[inc]=(ia,ih,oa,oh)
            if (ia or ih) and not (was[0] or was[1]):self.event(f'INNER self-contact detector inc={inc}, hitting={ih}; pair requires signed-gap validation')
            if (oa or oh) and not (was[2] or was[3]):self.event(f'OUTER self-contact detector inc={inc}, hitting={oh}; pair requires signed-gap validation')
            return
        m=RING.search(line)
        if m:
            side,inc,node,row,col,z=m.groups();key=(int(inc),side)
            if int(row) not in self.rings[key]:
                self.rings[key].add(int(row))
                self.event(f'{side} potential contact ring: inc={inc}, axial row={row}, distinct rows={len(self.rings[key])}; geometric confirmation still needed')
            return
        if 'NEWRAP CONVERGES' in line and self.cur is not None:
            latest=self.data[self.cur].get(self.iter)
            passed=latest is not None and latest<=self.tol
            self.conv[self.cur]=passed
            self.event(f'Increment {self.cur} solver-converged at NR iter {self.iter}, residual={latest:.7g} (tol={self.tol:.5g}); strict pass={passed}' if latest is not None else f'Increment {self.cur}: solver convergence flag (no residual)')
        elif 'N-R ITRATIONS DID NOT CONVERGE' in line:
            self.stopped=True;self.stopped_reason='Newton convergence failure'
            self.event(f'NON-CONVERGENCE at increment {self.cur} / Newton iteration {self.iter}.')
        elif 'Fortran runtime error' in line or 'Segmentation fault' in line:
            self.stopped=True; self.stopped_reason='Solver runtime error'; self.event('RUNTIME ERROR: '+line.strip())
        elif 'RUNNER TIME LIMIT' in line:
            self.stopped=True; self.stopped_reason='Execution-time limit reached (not Newton nonconvergence)'; self.event(line.strip())
    def svg(self, inc):
        d=self.data.get(inc,{})
        W,H=1050,500;l,r,t,b=95,35,47,75
        if not d:return ''
        xs=sorted(d); vals=[d[x] for x in xs]
        vv=[v for v in vals if v>0 and math.isfinite(v)]
        if not vv:return ''
        maxy=max(max(vv),self.tol)*1.5; miny=max(min(min(vv),self.tol)/2,1e-13)
        maxy=10**math.ceil(math.log10(maxy)); miny=10**math.floor(math.log10(miny))
        xmax=max(10,max(xs));xmin=0
        def X(v):return l+(W-l-r)*(v-xmin)/(xmax-xmin)
        def Y(v):return t+(H-t-b)*(math.log10(maxy)-math.log10(max(v,miny)))/(math.log10(maxy)-math.log10(miny))
        svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
             '<rect width="100%" height="100%" fill="#ffffff"/>',
             f'<text x="{l}" y="29" font-family="Arial" font-size="20" font-weight="bold" fill="#142c43">Increment {inc} — Newton–Raphson residual history</text>']
        for power in range(math.floor(math.log10(miny)),math.ceil(math.log10(maxy))+1):
            value=10**power;y=Y(value)
            svg.extend([f'<line x1="{l}" y1="{y:.1f}" x2="{W-r}" y2="{y:.1f}" stroke="#e4eaf2"/>',f'<text x="{l-12}" y="{y+5:.1f}" text-anchor="end" font-size="13" font-family="Arial" fill="#52616e">1e{power}</text>'])
        ticks=max(1,xmax//10)
        for it in range(0,xmax+1,ticks):
            x=X(it);svg.append(f'<line x1="{x:.1f}" y1="{t}" x2="{x:.1f}" y2="{H-b}" stroke="#f0f3f7"/><text x="{x:.1f}" y="{H-b+22}" text-anchor="middle" fill="#52616e" font-family="Arial" font-size="12">{it}</text>')
        yy=Y(self.tol)
        svg.append(f'<line x1="{l}" y1="{yy:.2f}" x2="{W-r}" y2="{yy:.2f}" stroke="#b3213d" stroke-width="2" stroke-dasharray="9,6"/><text x="{W-r-10}" y="{yy-9:.1f}" text-anchor="end" fill="#b3213d" font-family="Arial" font-size="13">Tolerance = {self.tol:.3g}</text>')
        pts=' '.join(f'{X(x):.2f},{Y(d[x]):.2f}' for x in xs if d[x]>0)
        svg.append(f'<polyline fill="none" stroke="#2865ad" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" points="{pts}"/>')
        for x in xs:
            v=d[x]
            if v<=0:continue
            svg.append(f'<circle cx="{X(x):.1f}" cy="{Y(v):.1f}" r="3.3" fill="#2865ad"><title>Iteration {x}: residual {v:.8g}</title></circle>')
        lx=xs[-1];lv=d[lx]
        svg.append(f'<text x="{l}" y="{H-15}" font-size="15" font-family="Arial" fill="#142c43">Newton–Raphson iteration</text>')
        svg.append(f'<text x="{W-r}" y="{H-15}" text-anchor="end" font-size="14" font-family="Arial" fill="#52616e">Last: iter {lx}, residual {lv:.6g} ({lv/self.tol:.1f}× tolerance)</text>')
        svg.append(f'<text transform="translate(23 {H/2}) rotate(-90)" text-anchor="middle" font-size="15" font-family="Arial" fill="#142c43">Unbalanced force ratio (log scale)</text>')
        svg.append('</svg>')
        return '\n'.join(svg)
    def render(self):
        if not self.data:return
        work=self.workdir;inc=self.cur if self.cur in self.data else max(self.data)
        for n in sorted(self.data):
            d=self.data[n]
            svg=self.svg(n)
            if svg:
                (work/f'increment_{n}_residual.svg').write_text(svg)
                with (work/f'increment_{n}_residual.csv').open('w',newline='') as fp:
                    wr=csv.writer(fp);wr.writerow(['newton_iteration','unbalanced_force_ratio','tolerance','strict_pass']);wr.writerows((it,val,self.tol,val<=self.tol) for it,val in sorted(d.items()))
        latest=self.data[inc];last=latest[max(latest)]
        inner,ih,outer,oh=self.contact.get(inc,(0,0,0,0))
        status=self.stopped_reason if self.stopped else ('Strictly converged' if self.conv.get(inc) else 'No verified convergence at this increment')
        summary=''.join(f'<tr><td>{n}</td><td>{max(d)}</td><td>{d[max(d)]:.7g}</td><td>{"PASS" if self.conv.get(n) else "Not verified"}</td></tr>' for n,d in sorted(self.data.items(),reverse=True)[:12])
        events=''.join('<li>'+html.escape(x)+'</li>' for x in reversed(self.events[-20:]))
        doc=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta http-equiv="refresh" content="8"><title>Cylinder FEM — Live residual monitor</title><style>
body{{font:16px Arial,sans-serif;background:#f5f7fb;color:#192a3c;margin:28px auto;max-width:1130px;padding:0 15px}}
h1{{font-size:27px}}.grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px}}.metric,.panel{{background:#fff;border:1px solid #e1e7ee;border-radius:9px;padding:15px;margin:12px 0}}.metric b{{font-size:24px;display:block;margin-top:5px}}.muted{{color:#697588;font-size:14px}}table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #edf0f4;padding:9px;text-align:left}}.panel svg{{width:100%;height:auto}}li{{margin:9px 0}}@media(max-width:700px){{.grid{{grid-template-columns:repeat(2,1fr)}}}}</style></head><body>
<h1>Cylinder FEM · residual and contact monitor</h1><p class="muted">Auto-refresh every 8 s. Last graph update: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}. This page monitors files only; it does not keep the executable running.</p>
<div class="grid"><div class="metric">Increment<b>{inc}</b></div><div class="metric">Newton iteration<b>{max(latest)}</b></div><div class="metric">Residual<b>{last:.4g}</b></div><div class="metric">Tolerance<b>{self.tol:.4g}</b></div></div>
<p><strong>Status:</strong> {html.escape(status)} · <strong>Last strict pass:</strong> {max((k for k,v in self.conv.items() if v), default='none')} · <strong>Inner:</strong> flag {inner}, hitting {ih} · <strong>Outer:</strong> flag {outer}, hitting {oh}</p>
<div class="panel">{self.svg(inc)}<p class="muted">Only Newton iteration residuals from the numerical log are plotted. The red dashed line represents the input convergence criterion.</p><p><a href="increment_{inc}_residual.svg">Open SVG plot</a> · <a href="increment_{inc}_residual.csv">Download CSV</a></p></div>
<div class="panel"><h2>Recent increments</h2><table><thead><tr><th>Increment</th><th>Last NR iteration</th><th>Residual</th><th>Convergence</th></tr></thead><tbody>{summary}</tbody></table></div>
<div class="panel"><h2>Contact / solver events</h2><ul>{events or '<li>No self-contact event recorded.</li>'}</ul><p class="muted">Additional axial rows are candidate rings, not validated contact reactions. A recorded detector event does not establish a physically admissible contact force.</p></div></body></html>'''
        tmp=work/'monitor.html.tmp';tmp.write_text(doc);os.replace(tmp,work/'monitor.html');self.last_render=time.monotonic();self.changed=False

def main():
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--follow',type=pathlib.Path);mode.add_argument('--replay',type=pathlib.Path)
    args=ap.parse_args();p=args.follow or args.replay; watcher=Watch(p.parent);offset=0;partial=''
    print(f'Monitoring {p} ; open {p.parent/"monitor.html"} in your browser.',flush=True)
    while True:
        try:
            with p.open('r',errors='replace') as fp:
                fp.seek(offset);chunk=fp.read();offset=fp.tell()
            dat=partial+chunk; lines=dat.splitlines(keepends=True);partial=''
            for line in lines:
                if line.endswith(('\n','\r')): watcher.process(line)
                else:partial=line
            if watcher.changed and (time.monotonic()-watcher.last_render>3 or args.replay):watcher.render()
            if args.replay:
                watcher.render();break
        except FileNotFoundError:
            if args.replay:sys.exit('Missing log file')
        time.sleep(2)

if __name__=='__main__': main()

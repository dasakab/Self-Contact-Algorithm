#!/usr/bin/env python3
"""Live, on-demand, dependency-free FEM dashboard for MSYS2 UCRT64.

Run from the cylinder project folder (with solver_live_windows.log):
    python3 cylinder_live_dashboard_windows.py --log solver_live_windows.log --port 8765

This serves its own HTML dashboard and reads the log on EACH browser request.
It never starts, stops, or changes the numerical solver; it only visualises logs.
"""
import argparse
import html
import math
import os
from pathlib import Path
import re
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from datetime import datetime
from urllib.parse import urlparse
import subprocess
import time

FLOAT = r'([-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?|NaN|Inf|Infinity|-Inf|-Infinity)'
STEP = re.compile(r'ADAPT_STEP_START\s+(\d+)\s+' + FLOAT + r'\s+' + FLOAT, re.I)
ITER = re.compile(r'INC NO\s*=\s*(\d+)\s+ITR NO\s*=\s*(\d+)', re.I)
RES = re.compile(r'UNBALANCE FORCE\s*=\s*' + FLOAT, re.I)
TOL = re.compile(r'TOLEC\s*=\s*' + FLOAT, re.I)
CONTACT = re.compile(r'SELF_CONTACT_STATUS inc inner outer:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', re.I)
RING = re.compile(r'SC_RING_(IN|OUT)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+'+FLOAT, re.I)
ACCEPT = re.compile(r'ADAPT_ACCEPT\s+(\d+)\s+' + FLOAT + r'\s+' + FLOAT + r'\s+(\d+)\s+(\d+)\s+(\d+)', re.I)
RETRY = re.compile(r'ADAPT_RETRY\s+(\d+)\s+', re.I)
ISSUE = re.compile(r'Fortran runtime error|Error termination|Segmentation fault|SIG(?:FPE|SEGV|ABRT)|FEM_NAN|N-R ITRATIONS DID NOT CONVERGE|malloc\(\):|Floating point exception', re.I)

def num(s):
    try:return float(s.replace('D','E').replace('d','e'))
    except ValueError:return float('nan')

def solver_status():
    # The check is informational only; it never affects the solver process.
    if os.name != 'nt': return None
    try:
        r=subprocess.run(['tasklist.exe', '/FI', 'IMAGENAME eq solver_windows.exe'], capture_output=True, text=True, timeout=3)
        if r.returncode == 0:
            return 'running' if re.search(r'(?im)^solver_windows\.exe\s+', r.stdout) else 'not running'
    except (OSError, subprocess.TimeoutExpired):pass
    return None

def load(log):
    if not log.exists():return None
    st=log.stat()
    # Limit read volume for long running calculations. Tail contains current history.
    with log.open('rb') as f:
        f.seek(max(0,st.st_size-16_000_000)); raw=f.read()
    lines=raw.decode('utf-8','replace').splitlines()
    current_inc=None; current_it=None
    values={}; last_tol=0.0005; latest_t=None; dt=None; contact={}; rows={}; accepts={}; errors=[]
    for ln in lines:
        if m:=STEP.search(ln):current_inc=int(m[1]);latest_t=num(m[2]);dt=num(m[3]);values.setdefault(current_inc,{})
        if m:=ITER.search(ln):current_inc=int(m[1]);current_it=int(m[2]);values.setdefault(current_inc,{})
        if m:=RETRY.search(ln):values[int(m[1])]={}
        if m:=RES.search(ln):
            if current_inc is not None and current_it is not None:values[current_inc][current_it]=num(m[1])
        if m:=TOL.search(ln):
            q=num(m[1]);
            if math.isfinite(q) and q>0:last_tol=q
        if m:=CONTACT.search(ln):contact[int(m[1])]=tuple(int(m[i]) for i in range(2,6))
        if m:=RING.search(ln):rows.setdefault((int(m[2]),m[1].upper()),set()).add(int(m[4]))
        if m:=ACCEPT.search(ln):accepts[int(m[1])]=(num(m[2]),num(m[3]),int(m[4]),int(m[5]),int(m[6]))
        if ISSUE.search(ln):errors.append(ln.strip()[:240])
    return dict(size=st.st_size,mtime=st.st_mtime,values=values,current=current_inc,it=current_it,
                time=latest_t,dt=dt,tol=last_tol,contacts=contact,rows=rows,accepts=accepts,
                errors=errors[-4:],tail=lines[-16:])

def curve(data,tol):
    pts=[(it,r) for it,r in sorted(data.items()) if math.isfinite(r) and r>0]
    if not pts:return '<p>No finite positive residuals to plot for this increment.</p>'
    W,H=920,360; L,R,T,B=85,20,30,58
    lower=max(min([tol]+[r for _,r in pts])/3,1e-14)
    upper=max([tol]+[r for _,r in pts])*1.5
    low=math.floor(math.log10(lower));high=math.ceil(math.log10(upper))
    if high<=low:high=low+1
    xmax=max(it for it,_ in pts)
    def xp(x):return L+(W-L-R)*(x-1)/max(1,xmax-1)
    def yp(v):return T+(H-T-B)*(high-math.log10(max(v,10**low)))/(high-low)
    out=[f'<svg role="img" aria-label="Newton residual chart" viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg">']
    for n in range(low,high+1):
        y=yp(10**n)
        out+=[f'<line x1="{L}" x2="{W-R}" y1="{y:.1f}" y2="{y:.1f}" stroke="#cbd5e1" opacity="0.6"/>',
              f'<text x="{L-8}" y="{y+5:.1f}" font-size="13" text-anchor="end" fill="#53637a">1e{n}</text>']
    y=yp(tol)
    out.append(f'<line x1="{L}" y1="{y:.1f}" x2="{W-R}" y2="{y:.1f}" stroke="#d43b54" stroke-width="2" stroke-dasharray="8 6"/>')
    out.append(f'<text x="{W-R-2}" y="{y-8:.1f}" text-anchor="end" font-size="12" fill="#c42543">Tolerance {tol:.4g}</text>')
    out.append('<polyline fill="none" stroke="#2563eb" stroke-width="2.5" points="'+' '.join(f'{xp(x):.1f},{yp(v):.1f}' for x,v in pts)+'"/>')
    out.append(f'<text x="{L}" y="{H-17}" font-size="14" fill="#334155">Newton iteration (1–{xmax})</text>')
    out.append('</svg>')
    return ''.join(out)

def render(log):
    data=load(log)
    now=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if data is None:
        body=f'<h1>Cylinder FEM — live residual monitor</h1><div class="warn">No log found: <code>{html.escape(str(log))}</code>. Start the numerical run or set <code>--log</code> to the right filename.</div>'
    else:
        age=max(0,int(time.time()-data['mtime']))
        process=solver_status()
        inc=data['current']; d=data['values'].get(inc,{}) if inc is not None else {}
        if not d and data['values']:
            inc=max(data['values']);d=data['values'][inc]
        latest=next(reversed(d.values())) if d else None
        latest_it=next(reversed(d)) if d else None
        nan=any(not math.isfinite(v) for v in d.values())
        msg=''
        if data['errors']:msg='Runtime / convergence message detected in log: '+data['errors'][-1]
        elif nan:msg='NaN/Inf was reported in the current residual history. This is not convergence.'
        elif process=='not running':msg='No solver_windows.exe process found. The solver appears stopped.'
        elif age>30:msg=f'No new log output for {age}s. The solver may be busy in factorisation, or stopped.'
        elif process=='running':msg='Solver process detected; this view reads fresh log data on every refresh.'
        else:msg='Monitoring the log. A fresh log line is needed to show numerical progress.'
        last_acc=max(data['accepts']) if data['accepts'] else None
        ct=data['contacts'].get(inc,(0,0,0,0)) if inc is not None else (0,0,0,0)
        rings_in=sorted(data['rows'].get((inc,'IN'),set())) if inc is not None else []
        rings_out=sorted(data['rows'].get((inc,'OUT'),set())) if inc is not None else []
        recent=''.join(f'<tr><td>{n}</td><td>{v[2]}</td><td>{v[0]:.9g}</td><td>{v[1]:.5g}</td><td>{v[3]}</td><td>{v[4]}</td></tr>' for n,v in sorted(data['accepts'].items(),reverse=True)[:12])
        details='\n'.join(data['tail'])
        body=f'''<h1>Cylinder FEM — live residual monitor</h1>
<div class="{ 'warn' if (nan or data['errors'] or process=='not running') else 'info'}">{html.escape(msg)}</div>
<div class="metrics"><div><small>Current increment</small><b>{inc if inc is not None else '—'}</b></div><div><small>Latest Newton iteration</small><b>{latest_it if latest_it is not None else '—'}</b></div><div><small>Latest residual</small><b>{str(round(latest,9)) if latest is not None else '—'}</b></div><div><small>Tolerance</small><b>{data['tol']:.3g}</b></div></div>
<p>Log last updated: <strong>{datetime.fromtimestamp(data['mtime']).strftime('%H:%M:%S')}</strong> ({age}s ago) · Size: {data['size']:,} bytes · Latest accepted increment: <strong>{last_acc if last_acc is not None else 'none in log'}</strong> · Time step: {data['dt'] if data['dt'] is not None else 'not in log'}</p>
<section><h2>Newton–Raphson residual vs iteration</h2>{curve(d,data['tol'])}</section>
<section><h2>Self-contact detectors</h2><p>Inner: flag={ct[0]}, hitting={ct[1]}, candidate axial rows={html.escape(str(rings_in))}.</p><p>Outer: flag={ct[2]}, hitting={ct[3]}, candidate axial rows={html.escape(str(rings_out))}.</p><p class="muted">Different axial rows are candidate rings, not independently verified physical contacts.</p></section>
<section><h2>Accepted increments</h2><table><tr><th>Increment</th><th>Newton iterations</th><th>Time</th><th>Δt</th><th>Inner flag</th><th>Outer flag</th></tr>{recent if recent else '<tr><td colspan="6">No ADAPT_ACCEPT records in current log</td></tr>'}</table></section>
<section><h2>Latest solver log lines</h2><pre>{html.escape(details)}</pre></section>'''
    return f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="5"><title>Cylinder FEM — live monitor</title><style>
body{{font:16px system-ui,Segoe UI,Arial;background:#f5f7fc;color:#1e293b;margin:25px auto;padding:0 20px;max-width:1020px}}h1{{font-size:25px}}h2{{font-size:19px}}section{{margin:16px 0;background:white;border:1px solid #dde4ec;border-radius:10px;padding:20px}}small{{display:block;color:#64748b;margin-bottom:8px}}b{{font-size:22px}}.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:16px 0}}.metrics div{{background:white;border:1px solid #dde4ec;border-radius:8px;padding:14px}}.warn,.info{{padding:14px;border-radius:8px;background:#fff4e4;border:1px solid #f3d09c}}.info{{background:#edf5ff;border-color:#b6d1f4}}.muted{{font-size:14px;color:#64748b}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;border-bottom:1px solid #e5e7eb;padding:7px}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;background:#f8fafc;padding:10px}}@media(max-width:700px){{.metrics{{grid-template-columns:repeat(2,1fr)}}}}
</style></head><body><p class="muted">On-demand log view · refreshed: {now} · automatic browser refresh every 5 seconds · no browser caching</p>{body}<p class="muted">This dashboard does not run the Fortran solver or fix convergence failures.</p></body></html>'''.encode('utf-8')

class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        path=urlparse(self.path).path
        if path not in ('/','/monitor.html','/live'):
            self.send_error(404,'Use /monitor.html');return
        response=render(self.server.log_path)
        self.send_response(200)
        self.send_header('Content-Type','text/html; charset=utf-8')
        self.send_header('Content-Length',str(len(response)))
        self.send_header('Cache-Control','no-store, no-cache, must-revalidate, max-age=0')
        self.send_header('Pragma','no-cache')
        self.send_header('Expires','0')
        self.end_headers();self.wfile.write(response)
    def log_message(self, fmt, *args):
        pass

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--log',default='solver_live_windows.log',help='Numerical solver output logfile')
    ap.add_argument('--port',type=int,default=8765)
    args=ap.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    server.log_path=Path(args.log).resolve()
    print('Live dashboard: http://localhost:%d/monitor.html' % args.port,flush=True)
    print('Reading log:',server.log_path,flush=True)
    print('Refreshes every 5 seconds; press Ctrl+C to stop server.',flush=True)
    try:server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:pass
    finally:server.server_close()

if __name__=='__main__':main()

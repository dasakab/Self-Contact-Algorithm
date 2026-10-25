#!/usr/bin/env python3
"""Display measured Newton residuals, adaptive dt decisions and separate contacts.
Follow a locally running solver with: python3 adaptive_dashboard.py --follow adaptive_run.log
Replay a completed log with:    python3 adaptive_dashboard.py --replay adaptive_run.log
The HTML automatically refreshes; it cannot launch the Fortran solver.
"""
import argparse, re, time, pathlib, html, collections, csv, os
from residual_dashboard import Watch
F=r'([-+\d.eEdD]+)'
step=re.compile(r'ADAPT_STEP_START\s+(\d+)\s+'+F+r'\s+'+F+r'\s+(\d+)')
accept=re.compile(r'ADAPT_ACCEPT\s+(\d+)\s+'+F+r'\s+'+F+r'\s+(\d+)\s+(\d+)\s+(\d+)')
nextstep=re.compile(r'ADAPT_NEXT\s+(\d+)\s+'+F+r'\s+(\d+)')
retry=re.compile(r'ADAPT_RETRY\s+(\d+)\s+(\d+)\s+'+F+r'\s+'+F)
cut=re.compile(r'ADAPT_CUT_APPLIED\s+(\d+)\s+'+F+r'\s+(\d+)')
from residual_dashboard import f
class AdaptiveWatch(Watch):
    def __init__(self,workdir):
        super().__init__(workdir)
        self.timeline=[]; self.attempts=collections.defaultdict(int); self.active_dt={}; self.accepted={}; self.saved_attempts=[]
    def process(self,line):
        if m:=step.search(line):
            n,t,dt,streak=m.groups();n=int(n);self.stopped=False;self.stopped_reason='';self.active_dt[n]=float(dt);self.timeline.append((n,'STEP START',f't={t}  dt={dt}  good streak={streak}'));self.changed=True
            self.event(f'Increment {n} start: t={t}, Δt={dt}')
        elif m:=retry.search(line):
            n,it,dt,r=m.groups();n=int(n); self.timeline.append((n,'RETRY',f'iteration={it}, dt={dt}, residual={r}'));self.event(f'Increment {n}: retry required at iter {it}; residual {r}')
            if n in self.data and self.data[n]:
                old=dict(self.data[n]);number=self.attempts[n]+1
                with (self.workdir/f'increment_{n}_attempt_{number}_residual.csv').open('w',newline='') as fh:
                    cw=csv.writer(fh);cw.writerow(['iteration','residual']);cw.writerows(sorted(old.items()))
                svg=self.svg(n)
                (self.workdir/f'increment_{n}_attempt_{number}.svg').write_text(svg)
                self.data[n].clear()
            self.attempts[n]+=1
        elif m:=cut.search(line):
            n,dt,num=m.groups();self.active_dt[int(n)]=float(dt);self.timeline.append((int(n),'TIME CUT',f'new dt={dt}, cut #{num}'));self.changed=True
        elif m:=accept.search(line):
            n,t,dt,it,ih,oh=m.groups();n=int(n);self.accepted[n]=(float(t),float(dt),int(it),int(ih),int(oh))
            self.timeline.append((n,'ACCEPTED',f't={t}, Δt={dt}, NR={it}, inner={ih}, outer={oh}'));self.changed=True
        elif m:=nextstep.search(line):
            n,dt,streak=m.groups();self.timeline.append((int(n),'NEXT Δt',f'planned dt={dt} good streak={streak}'));self.changed=True
        if 'malloc():' in line or 'SIGABRT' in line:
            self.stopped=True
            self.stopped_reason='Heap corruption / process abort (not Newton nonconvergence)'
            self.event(f'RUNTIME MEMORY FAILURE: {line.strip()}')
        super().process(line)
    def render(self):
        if not self.data:return
        super().render()
        p=self.workdir/'monitor.html'
        s=p.read_text()
        latest=max(self.data)
        dt=self.active_dt.get(latest)
        rows=''.join(f'<tr><td>{n}</td><td>{html.escape(st)}</td><td>{html.escape(detail)}</td></tr>' for n,st,detail in reversed(self.timeline[-35:]))
        latestpass=max((n for n,ok in self.conv.items() if ok),default=None)
        cp='not available'
        pp=self.workdir/'RESTART'/'PARA'
        if pp.exists():
            try:
                raw=pp.read_text().split();cp=f'increment {int(raw[0])}, t={float(raw[1]):.9g}'
            except (ValueError,IndexError):pass
        extra=f'''<div class="panel"><h2>Adaptive time-step controller</h2><p><b>Current Δt:</b> {dt if dt is not None else 'Not found in current log'} · <b>Saved checkpoint:</b> {cp} · <b>Latest strict pass in log:</b> {latestpass}</p><p class="muted">Increase ×1.25 after two fast convergences (≤60 Newton iterations). Reduce ×0.8 at ≥80, ×0.5 at ≥100. Retry failures with half Δt down to 1e−7. Independent inner/outer detectors are tracked above. Active self-contact limits future growth to Δt ≤1e−5.</p><table><thead><tr><th>Increment</th><th>Decision</th><th>Details</th></tr></thead><tbody>{rows}</tbody></table></div>'''
        p.write_text(s.replace('</body>',extra+'</body>'))
        # The parent render stores a monotonic timestamp for incremental refresh.
def feed(w,source):
    try:
        with open(source,encoding='utf-8',errors='replace') as fh:
            for line in fh:w.process(line)
    except FileNotFoundError:pass

def main():
    ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--follow',type=pathlib.Path);g.add_argument('--replay',type=pathlib.Path)
    ap.add_argument('--history',type=pathlib.Path,default=None)
    args=ap.parse_args();log=args.follow or args.replay
    w=AdaptiveWatch(log.parent)
    if args.history:feed(w,args.history)
    if args.replay:
        feed(w,log);w.render(); print('Dashboard: '+str(log.parent/'monitor.html'));return
    off=0; partial='';last=0
    while True:
        try:
            with log.open(errors='replace') as fh:
                fh.seek(off); data=fh.read();off=fh.tell()
            lines=(partial+data).splitlines(keepends=True);partial=''
            for line in lines:
                if line.endswith(('\n','\r')): w.process(line)
                else:partial=line
            if w.changed and time.monotonic()-last>5:
                w.render();last=time.monotonic()
        except FileNotFoundError:pass
        time.sleep(2)
if __name__=='__main__':main()

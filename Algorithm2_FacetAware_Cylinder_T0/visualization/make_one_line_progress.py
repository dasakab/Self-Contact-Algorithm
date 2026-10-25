#!/usr/bin/env python3
"""Print one-sentence status per accepted FEM increment, plus pending attempt."""
import re,sys
from pathlib import Path
p=Path(sys.argv[1]) if len(sys.argv)>1 else Path('trial.log')
lines=p.read_text(errors='replace').splitlines()
contacts={}; rows={}; rs={}; vals={}; accepts=[]
cur=0
for line in lines:
    m=re.search(r'UNBALANCE FORCE\s*=\s*([\d.eE+-]+)',line)
    if m and cur:
        try: vals[cur]=float(m.group(1).replace('D','E'))
        except: pass
    m=re.search(r'INC NO\s*=\s*(\d+)\s+ITR NO\s*=\s*(\d+)',line)
    if m: cur=int(m.group(1))
    m=re.search(r'SELF_CONTACT_STATUS inc inner outer:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',line)
    if m: contacts[int(m.group(1))]=(int(m.group(3)),int(m.group(5)))
    m=re.search(r'SC_RING_(IN|OUT)\s+(\d+)\s+\d+\s+(\d+)',line)
    if m:rows.setdefault(int(m.group(2)),{'IN':set(),'OUT':set()})[m.group(1)].add(int(m.group(3)))
    m=re.search(r'ADAPT_ACCEPT\s+(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+(\d+)',line)
    if m:
        inc=int(m.group(1));t=float(m.group(2));dt=float(m.group(3));it=int(m.group(4));
        rr=vals.get(inc)
        ins,outs=contacts.get(inc,(None,None))
        row=rows.get(inc,{'IN':set(),'OUT':set()})
        desc=f'Increment {inc}: converged in {it} Newton iterations at t={t:.9f}, dt={dt:.3g}, residual={rr:.8g}' if rr is not None else f'Increment {inc}: accepted after {it} iterations at t={t:.9f}, dt={dt:.3g}'
        desc+=f'; inner hits={ins}, rows={sorted(row["IN"])}, outer hits={outs}, rows={sorted(row["OUT"])}.'
        accepts.append(desc)
for s in accepts:print(s)
if cur and (not accepts or cur>int(re.search(r'Increment (\d+)',accepts[-1]).group(1))):
    print(f'Increment {cur}: attempted, NOT accepted in available log; latest residual={vals.get(cur,"not recorded")}.')

#!/usr/bin/env python3
"""Follow solver log. Record contact rows for inner and outer surfaces independently.
An event is accepted only when the corresponding increment is accepted by the FEM solver.
Detected candidates still require physical verification (normal gap / non-adjacent facets).
"""
import argparse,collections,csv,time,re,pathlib
p=argparse.ArgumentParser();p.add_argument('log',nargs='?',default='solver_live.log');p.add_argument('--follow',action='store_true');args=p.parse_args()
log=pathlib.Path(args.log); folder=log.parent; rows=collections.defaultdict(lambda:{'IN':collections.defaultdict(set),'OUT':collections.defaultdict(set)})
first={}; emitted=set(); accepted={}; counts={}
pat_ring=re.compile(r'SC_RING_(IN|OUT)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
pat_counts=re.compile(r'SELF_CONTACT_STATUS inc inner outer:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
pat_accept=re.compile(r'ADAPT_ACCEPT\s+(\d+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)\s+(\d+)')
out=folder/'contact_ring_events.csv'
with out.open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['increment','time','status','surface','axial_row','node_count_in_row','inner_hits','outer_hits'])
    def emit(n,accepted_now=False):
        for surface in ('IN','OUT'):
            data=rows[n][surface]
            if not data:continue
            for row,nodes in sorted(data.items()):
                key=(n,surface,row)
                if key in emitted:continue
                emitted.add(key)
                if surface not in first:first[surface]=row
                phase=('NEW_AXIAL_ROW_CANDIDATE' if row!=first[surface] else 'FIRST_ROW_OR_SAME_ROW')
                if surface=='OUT':phase='OUTER_'+phase
                timev=accepted.get(n,'') if accepted_now else ''
                hi,ho=counts.get(n,(None,None))
                w.writerow([n,timev,phase,surface,row,len(nodes),hi,ho]);f.flush()
                print(f'INC {n}, {surface} contact, row {row}, {len(nodes)} nodes: {phase}',flush=True)
                if 'NEW_AXIAL_ROW_CANDIDATE' in phase and accepted_now:
                    (folder/'SECOND_RING_CANDIDATE.txt').write_text(f'Increment {n}: {surface} new axial ring candidate at row {row}; solver accepted at time {timev}. Verify signed gaps and contact forces.\n')
    offset=0;part=''
    while True:
        try:
            with log.open(errors='replace') as fd:
                fd.seek(offset);fresh=fd.read();offset=fd.tell()
            lines=(part+fresh).splitlines(keepends=True);part=''
            for line in lines:
                if not line.endswith(('\n','\r')):part=line;continue
                m=pat_ring.search(line)
                if m:
                    surf,n,node,row,_=m.groups();rows[int(n)][surf][int(row)].add(int(node))
                m=pat_counts.search(line)
                if m:
                    n,ia,ih,oa,oh=map(int,m.groups());counts[n]=(ih,oh)
                m=pat_accept.search(line)
                if m:
                    n=int(m[1]);accepted[n]=m[2];emit(n,True)
        except FileNotFoundError:pass
        if not args.follow:break
        time.sleep(3)

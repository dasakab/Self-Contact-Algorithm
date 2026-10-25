#!/usr/bin/env python3
"""Configure numerical input without changing the Fortran physics.

This script edits exactly the 6th record of UOPT1/OPT_IEXEC and optionally
the first 3 records of UOPT1/O_OEXEC. Use through launcher scripts.
"""
import argparse
from pathlib import Path
import math

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
ap.add_argument('--mode',choices=['zero','restart'],required=True)
ap.add_argument('--end-time',type=float,default=None)
ap.add_argument('--initial-dt',type=float,default=None)
ap.add_argument('--max-iterations',type=int,default=None)
args=ap.parse_args()
root=Path(args.root).resolve()
file=root/'UOPT1/OPT_IEXEC'
fields=file.read_text().splitlines()
if len(fields)<6:raise SystemExit('OPT_IEXEC must have 6 records')
fields[5]=f'{0 if args.mode=="zero" else 1:12d}'
file.write_text('\n'.join(fields)+'\n')
cfg=root/'UOPT1/O_OEXEC'
vals=cfg.read_text().splitlines()
if len(vals)<8:raise SystemExit('O_OEXEC must have at least 8 records')
if args.end_time is not None:
    if not (math.isfinite(args.end_time) and args.end_time>0):raise SystemExit('end time must be positive and finite')
    vals[0]=f'{args.end_time:24.15E}'
if args.max_iterations is not None:
    if args.max_iterations<1:raise SystemExit('iteration maximum must be positive')
    vals[1]=f'{args.max_iterations:12d}'
if args.initial_dt is not None:
    if not (math.isfinite(args.initial_dt) and args.initial_dt>0):raise SystemExit('time step must be positive and finite')
    vals[2]=f'{args.initial_dt:24.15E}'
cfg.write_text('\n'.join(vals)+'\n')
print(f'Input configured: mode={args.mode} (restart flag={fields[5].strip()}); end_time={vals[0].strip()}, max_iterations={vals[1].strip()}, starting dt={vals[2].strip()}')

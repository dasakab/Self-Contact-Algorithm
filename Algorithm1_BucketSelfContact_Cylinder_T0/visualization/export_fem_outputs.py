#!/usr/bin/env python3
"""Export this Fortran FEM solver's OUT + UBIN1 + RESTART results to ParaView.

No third-party modules needed for VTU/PVD conversion.
Optional --gif needs matplotlib and Pillow (available in MSYS2 UCRT64).

Usage:
    python export_fem_outputs.py --root .
    python export_fem_outputs.py --root . --gif

Reads the complete per-increment blocks in OUT/*AFTER_DEFORM. Optionally adds
an explicitly named checkpoint with --checkpoint, or undeformed reference
configuration with --include-initial. Missing increments are NEVER interpolated.
"""
import argparse
import math
from pathlib import Path
import re
import sys
from collections import defaultdict
from xml.sax.saxutils import escape

INC_RE = re.compile(r'^\s*INCREMENT\s+NUMBER:\s*(\d+)', re.IGNORECASE)
STEP_RE = re.compile(r'ADAPT_STEP_START\s+(\d+)\s+([\d.EeDd+-]+)\s+([\d.EeDd+-]+)')


def number(value):
    return float(value.replace('D','E').replace('d','e'))


def read_indexed_file(filename, ncomp):
    output = {}
    for ln in filename.read_text(encoding='utf-8',errors='replace').splitlines():
        parts = ln.split()
        if len(parts) < ncomp+1:
            continue
        try:
            i=int(parts[0]); values=tuple(number(x) for x in parts[1:1+ncomp]);
            if i >= 1 and all(math.isfinite(x) for x in values):
                output[i]=values
        except (ValueError, OverflowError):
            continue
    return output


def read_increment_blocks(filename, ncomp):
    if not filename.is_file():
        print('Missing:', filename)
        return {}
    output={}; current=None
    for ln in filename.read_text(encoding='utf-8', errors='replace').splitlines():
        match=INC_RE.search(ln)
        if match:
            current=int(match.group(1)); output[current]={}
            continue
        if current is None:
            continue
        parts=ln.split()
        if len(parts) < ncomp+1:
            continue
        try:
            i=int(parts[0]); vals=tuple(number(x) for x in parts[1:ncomp+1]);
            if i >= 1 and all(math.isfinite(v) for v in vals):
                output[current][i]=vals
        except (ValueError,OverflowError):
            continue
    return output


def read_connectivity(filename):
    output={}
    for ln in filename.read_text(errors='replace').splitlines():
        parts=ln.split()
        if len(parts)<9: continue
        try:
            idx=int(parts[0]); nodes=tuple(int(x) for x in parts[1:9])
            if idx>=1 and all(x>=1 for x in nodes):output[idx]=nodes
        except ValueError: continue
    return output


def read_time_table(root):
    times={}
    for filename in [root/'continuation.log',root/'solver_live.log',root/'solver_live_windows.log',
                     root/'diagnostic_run.log',root/'run.log',root/'OUT/EXEC']:
        if not filename.is_file(): continue
        for ln in filename.read_text(errors='replace').splitlines():
            m=STEP_RE.search(ln)
            if m:
                inc=int(m.group(1)); start=number(m.group(2)); dt=number(m.group(3))
                times[inc]=start+dt
    # Latest checkpoint time is authoritative for last accepted increment
    para=root/'RESTART/PARA'
    if para.is_file():
        p=para.read_text(errors='replace').split()
        try:
            times[int(p[0])]=number(p[1])
        except (ValueError, IndexError): pass
    return times


def write_vtu(filename, coords, conn, displacements, stresses, plastic, inc, time):
    node_ids=sorted(coords)
    elem_ids=sorted(conn)
    id_to_zero={nid:j for j,nid in enumerate(node_ids)}
    cells=[]; surviving_eids=[]
    for eid in elem_ids:
        if all(n in id_to_zero for n in conn[eid]):
            cells.append(tuple(id_to_zero[n] for n in conn[eid])); surviving_eids.append(eid)
    if not cells:raise ValueError('No valid elements mapped to coordinates')
    def fmt(x):return format(x,'.11g') if math.isfinite(x) else 'nan'
    def p_data(name,comps,writer):
        writer.write(f'        <DataArray type="Float64" Name="{name}" NumberOfComponents="{comps}" format="ascii">\n')
    with filename.open('w',encoding='utf-8',newline='\n') as w:
        w.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        w.write('  <UnstructuredGrid>\n')
        w.write(f'    <Piece NumberOfPoints="{len(node_ids)}" NumberOfCells="{len(cells)}">\n')
        w.write('      <Points>\n')
        p_data('Points',3,w)
        for nid in node_ids: w.write('          '+' '.join(map(fmt,coords[nid]))+'\n')
        w.write('        </DataArray>\n      </Points>\n      <Cells>\n')
        w.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        for c in cells: w.write('          '+' '.join(map(str,c))+'\n')
        w.write('        </DataArray>\n        <DataArray type="Int64" Name="offsets" format="ascii">\n          ')
        w.write(' '.join(str(8*(i+1)) for i in range(len(cells))))
        w.write('\n        </DataArray>\n        <DataArray type="UInt8" Name="types" format="ascii">\n          ')
        w.write(' '.join('12' for _ in cells)) # VTK_HEXAHEDRON
        w.write('\n        </DataArray>\n      </Cells>\n      <PointData Vectors="displacement_vector">\n')
        p_data('displacement_vector',3,w)
        for nid in node_ids:
            vals=displacements.get(nid,(float('nan'),)*3)
            w.write('          '+' '.join(map(fmt,vals))+'\n')
        w.write('        </DataArray>\n')
        p_data('displacement_magnitude',1,w)
        for nid in node_ids:
            d=displacements.get(nid,(float('nan'),)*3)
            w.write('          '+fmt(math.sqrt(sum(x*x for x in d)))+'\n')
        w.write('        </DataArray>\n')
        w.write('        <DataArray type="Int32" Name="original_node_id" format="ascii">\n          ')
        w.write(' '.join(map(str,node_ids)))
        w.write('\n        </DataArray>\n      </PointData>\n      <CellData Scalars="equivalent_stress_Pa">\n')
        for name,source in [('equivalent_stress_Pa',stresses),('equivalent_plastic_strain',plastic)]:
            p_data(name,1,w)
            for eid in surviving_eids:
                val=source.get(eid,(float('nan'),))[0]
                w.write('          '+fmt(val)+'\n')
            w.write('        </DataArray>\n')
        w.write('        <DataArray type="Int32" Name="original_element_id" format="ascii">\n          ')
        w.write(' '.join(map(str,surviving_eids)))
        w.write('\n        </DataArray>\n      </CellData>\n')
        w.write(f'      <FieldData>\n        <DataArray type="Int32" Name="increment" NumberOfTuples="1" format="ascii">{inc}</DataArray>\n')
        if time is not None:
            w.write(f'        <DataArray type="Float64" Name="time_s" NumberOfTuples="1" format="ascii">{fmt(time)}</DataArray>\n')
        w.write('      </FieldData>\n')
        w.write('    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n')
    return node_ids,surviving_eids


def get_slice_faces(initial,connectivity):
    face_sets=((0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    out=[]
    for eid,elem in connectivity.items():
        for fs in face_sets:
            q=tuple(elem[i] for i in fs)
            if all(n in initial and abs(initial[n][1])<1e-8 for n in q):
                out.append((eid,q));break
    return out


def create_gif(result_states, initial, connectivity, output):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        from PIL import Image
    except ImportError as e:
        print('GIF skipped: install matplotlib and Pillow. Detail:',e)
        return False
    slice_faces=get_slice_faces(initial, connectivity)
    if not slice_faces:
        print('GIF skipped: no mesh faces lie on reference y=0 plane')
        return False
    # constant framing, same plane and colours for all frames
    xmn=ymn=float('inf');xmx=ymx=float('-inf')
    for inc,tm,coords,_,_,_ in result_states:
        for _,quad in slice_faces:
            for nid in quad:
                x,z=coords[nid][0],coords[nid][2]
                xmn=min(xmn,x);xmx=max(xmx,x);ymn=min(ymn,z);ymx=max(ymx,z)
    xmn-=.1;xmx+=.1;ymn-=.1;ymx+=.1
    image_list=[]
    frames_folder=output.parent/'PNG_frames'
    frames_folder.mkdir(exist_ok=True)
    for inc,tm,coords,displacements,stresses,plastic in result_states:
        polys=[];v=[]
        for eid,quad in slice_faces:
            polys.append([(coords[n][0],coords[n][2]) for n in quad]);v.append(stresses.get(eid,(float('nan'),))[0]/1e6)
        fig,ax=plt.subplots(figsize=(6,10),facecolor='#596078')
        ax.set_facecolor('#596078')
        pc=PolyCollection(polys,array=__import__('numpy').asarray(v,dtype=float),cmap='coolwarm',edgecolors='none')
        pc.set_clim(0,1100)
        ax.add_collection(pc)
        mesh=PolyCollection(polys,facecolors='none',edgecolors=(1,1,1,.17),linewidths=.35)
        ax.add_collection(mesh)
        ax.set_xlim(xmn,xmx);ax.set_ylim(ymn,ymx)
        ax.set_aspect('equal');ax.axis('off')
        fig.colorbar(pc,ax=ax,shrink=.76,label='Equivalent stress (MPa)')
        title=f'Increment {inc}'
        if tm is not None:title+=f'  |  t={tm:.8f} s'
        fig.suptitle(title,color='white',fontsize=13,y=.98)
        fig.subplots_adjust(left=.07,right=.87,top=.94,bottom=.04)
        fn=frames_folder/f'frame_{inc:04d}.png'
        fig.savefig(fn,dpi=120,facecolor=fig.get_facecolor())
        plt.close(fig)
        with Image.open(fn) as im:
            image_list.append(im.convert('RGB').copy())
    if image_list:
        image_list[0].save(output,save_all=True,append_images=image_list[1:],duration=[1100]+[700]*(len(image_list)-2)+([1400] if len(image_list)>1 else []),loop=0,optimize=True)
        print('GIF:',output,'(',len(image_list),'actual saved frames)')
        return True
    return False


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root',default='.',help='folder holding UBIN1, OUT, RESTART')
    ap.add_argument('--gif',action='store_true',help='make cross-section GIF (requires matplotlib, Pillow)')
    ap.add_argument('--output',default='PARAVIEW_EXPORT',help='output folder name')
    ap.add_argument('--checkpoint',default=None,
                    help='Optional saved checkpoint folder, e.g. examples/checkpoints/INC35_VERIFIED. Never loaded by default.')
    ap.add_argument('--include-initial',action='store_true',
                    help='Export reference undeformed geometry at time zero with zero displacement/stress/strain; not a computed result.')
    args=ap.parse_args()
    root=Path(args.root).expanduser().resolve()
    out=root/args.output;out.mkdir(exist_ok=True)
    conn_file=root/'UBIN1/ELE_NODE_A';coord_file=root/'UBIN1/COORD_A'
    if not conn_file.is_file() or not coord_file.is_file():
        raise SystemExit('Expected UBIN1/ELE_NODE_A and UBIN1/COORD_A under '+str(root))
    conn=read_connectivity(conn_file)
    initial=read_indexed_file(coord_file,3)
    coords=read_increment_blocks(root/'OUT/COORD_A_AFTER_DEFORM',3)
    disp=read_increment_blocks(root/'OUT/DISP_A_AFTER_DEFORM',3)
    stress=read_increment_blocks(root/'OUT/CSTRESS_A_AFTER_DEFORM',1)
    plast=read_increment_blocks(root/'OUT/EQUI_PS_A_AFTER_DEFORM',1)
    times=read_time_table(root)
    states=[]
    if args.include_initial:
        # MREAD adds the body-A origin offset from UOPT1/MESH to UBIN1/COORD_A.
        # Apply that same rigid translation to the initial/reference VTU.
        mesh_options=(root/'UOPT1/MESH').read_text(errors='replace').splitlines()
        origin_vals=[]
        for line in reversed(mesh_options):
            try:
                v=[number(t) for t in line.split()]
                if len(v)==3:origin_vals=v;break
            except ValueError:pass
        if len(origin_vals)!=3:
            raise SystemExit('Cannot determine body-A placement offset in UOPT1/MESH for --include-initial')
        initial_placed={i:tuple(initial[i][k]+origin_vals[k] for k in range(3)) for i in initial}
        zeros_n={i:(0.0,0.0,0.0) for i in initial}
        zeros_e={i:(0.0,) for i in conn}
        states.append((0,0.0,initial_placed,zeros_n,zeros_e,zeros_e))
    ck=(root/args.checkpoint).resolve() if args.checkpoint else None
    if ck is not None and (ck/'COORD_A_DEFORM_CUR_INC').is_file():
        c=read_indexed_file(ck/'COORD_A_DEFORM_CUR_INC',3)
        d=read_indexed_file(ck/'DISP_A_DEFORM_CUR_INC',3)
        s=read_indexed_file(ck/'EQUI_STRESS_A_CUR_INC',1)
        p=read_indexed_file(ck/'EQUI_P_STRAIN_A_CUR_INC',1)
        inc=35;tm=None
        para=ck/'PARA'
        if para.is_file():
            try: inc0,tm0=para.read_text().split()[:2];inc=int(inc0);tm=number(tm0)
            except (ValueError,IndexError):pass
        if len(c)==len(initial) and len(d)==len(initial) and len(s)==len(conn) and len(p)==len(conn):
            states.append((inc,tm,c,d,s,p))
        else:
            print('Incomplete explicit checkpoint not included:',ck)
    if args.checkpoint and (ck is None or not (ck/'COORD_A_DEFORM_CUR_INC').exists()):
        raise SystemExit('Checkpoint not found: '+str(ck))
    for inc in sorted(coords):
        c=coords[inc];d=disp.get(inc,{});s=stress.get(inc,{});p=plast.get(inc,{})
        if len(c)!=len(initial) or len(d)!=len(initial) or len(s)!=len(conn) or len(p)!=len(conn):
            print(f'Skipping incomplete increment {inc}: nodes={len(c)}/{len(initial)}, disps={len(d)}, stress={len(s)}/{len(conn)}, epsp={len(p)}')
            continue
        states.append((inc,times.get(inc),c,d,s,p))
    # Keep latest unique inc
    states=list({s[0]:s for s in states}.values());states.sort(key=lambda s:s[0])
    if not states:raise SystemExit('No complete saved increments found. Check OUT files and checkpoint.')
    manifest=[]
    for inc,tm,c,d,s,p in states:
        fn=out/f'cylinder_inc{inc:04d}.vtu'
        write_vtu(fn,c,conn,d,s,p,inc,tm)
        peak=max(v[0] for v in s.values())/1e6
        print(f'Increment {inc}: {fn.name}, t={tm}, max equivalent stress={peak:.1f} MPa')
        manifest.append((inc,tm,fn.name))
    time_type='physical time in seconds' if all(t is not None for _,t,_ in manifest) else 'mixed: increment index when time unknown'
    pvd=out/'cylinder_timesteps.pvd'
    with pvd.open('w',encoding='utf-8') as w:
        w.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n')
        for inc,tm,filename in manifest:
            w.write(f'    <DataSet timestep="{tm if tm is not None else inc}" group="" part="0" file="{filename}"/>\n')
        w.write('  </Collection>\n</VTKFile>\n')
    with (out/'README.txt').open('w') as w:
        w.write('ParaView FEM time series export\n')
        w.write(f'Input: {root}\nSaved frames: {len(states)}\nPVD time axis: {time_type}\n')
        w.write('All values come from stored solver output, not interpolated.\n')
        w.write('To view: open cylinder_timesteps.pvd in ParaView, Apply, Slice with Normal=(0,1,0), use Camera +Z up looking along Y, choose equivalent_stress_Pa or displacement_magnitude.\n')
        w.write('Stress data in Pa, displacement/coordinates in model units (check solver input units).\n')
    print('ParaView time series:',pvd)
    if args.gif:create_gif(states,initial,conn,out/'cylinder_cross_section.gif')
    else:print('Add --gif to create a cross-section GIF with matplotlib + Pillow')
    actual_count = len(states) - (1 if args.include_initial and any(item[0] == 0 for item in states) else 0)
    print(f'Finished. Exported {actual_count} complete saved increments' +
          (' plus the initial undeformed reference' if args.include_initial else '') +
          '. No intermediate increments were invented.')


if __name__=='__main__':main()

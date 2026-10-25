#!/usr/bin/env python3
"""Detect nonincident triangle-triangle crossings in saved cylinder boundary mesh.
No force or contact validity is inferred from visual overlap alone.
"""
from pathlib import Path
import sys, vtk, numpy as np
from scipy.spatial import cKDTree
from vtk.util.numpy_support import vtk_to_numpy

fn=Path(sys.argv[1]); reader=vtk.vtkXMLUnstructuredGridReader();reader.SetFileName(str(fn));reader.Update(); grid=reader.GetOutput()
surf=vtk.vtkDataSetSurfaceFilter();surf.SetInputData(grid);surf.Update()
tri=vtk.vtkTriangleFilter();tri.SetInputData(surf.GetOutput());tri.Update();pol=tri.GetOutput()
pts=vtk_to_numpy(pol.GetPoints().GetData())
triangles=np.array([[pol.GetCell(i).GetPointId(j) for j in range(3)] for i in range(pol.GetNumberOfCells())],dtype=np.int32)
xyz=pts[triangles]
center=xyz.mean(axis=1)
rad=np.sqrt(((xyz-center[:,None])**2).sum(axis=2)).max(axis=1)
searchradius=rad.max()*2
if not np.isfinite(searchradius):raise RuntimeError('nonfinite mesh')
tree=cKDTree(center)
count=0; found=[]; checks=0
for i in range(len(triangles)):
    for j in tree.query_ball_point(center[i], rad[i]+rad.max()+1e-9):
        if j<=i:continue
        if np.linalg.norm(center[j]-center[i])>rad[i]+rad[j]+1e-10:continue
        if np.intersect1d(triangles[i],triangles[j]).size:continue
        checks+=1
        if vtk.vtkTriangle.TrianglesIntersect(*[xyz[i,k].tolist() for k in range(3)],*[xyz[j,k].tolist() for k in range(3)]):
            count+=1
            if len(found)<25:found.append((int(i),int(j),xyz[i].mean(axis=0).tolist(),xyz[j].mean(axis=0).tolist()))
print('SOURCE',fn.name,'TRIANGLES',len(triangles),'NONINCIDENT_TESTS',checks,'TRIANGLE_CROSSINGS',count,flush=True)
for i,j,a,b in found:print('CROSS',i,j,'CENTER1',a,'CENTER2',b)

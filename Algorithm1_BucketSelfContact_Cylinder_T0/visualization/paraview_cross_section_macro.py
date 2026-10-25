"""ParaView GUI Python Shell macro; open cylinder_timesteps.pvd first and Apply.
Execute via View > Python Shell > Run Script, or add to Tools > Manage Macros.
Keeps camera on X-Z cross-section, looks along Y, with Z as vertical.
"""
from paraview.simple import (GetActiveSource, Slice, Show, GetActiveViewOrCreate,
                             ColorBy, Hide, Render, ResetCamera)
source=GetActiveSource()
if source is None:
    raise RuntimeError('Open PARAVIEW_EXPORT/cylinder_timesteps.pvd and Apply first')
view=GetActiveViewOrCreate('RenderView')
section=Slice(registrationName='Cross section (y=0)',Input=source)
section.SliceType='Plane'
section.SliceType.Origin=[0.0,0.0,22.0]
section.SliceType.Normal=[0.0,1.0,0.0]
Hide(source,view)
disp=Show(section,view)
disp.Representation='Surface With Edges'
ColorBy(disp,('CELLS','equivalent_stress_Pa'))
view.CameraFocalPoint=[0.0,0.0,22.0]
view.CameraPosition=[0.0,-10.0,22.0]
view.CameraViewUp=[0.0,0.0,1.0]
view.CameraParallelProjection=1
ResetCamera(view)
Render()

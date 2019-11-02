
from paraview.simple import *


n=1009

s =  "%04d" % n
file_psi  = ["eps1e-2_4096_zoom_VTK/psi." + s + ".vtk"]
file_ener = ["eps1e-2_4096_zoom_VTK/ener." + s + ".vtk"]


#-- parameters --

hscale = 30.0  # scale fot the height of psi

emin = -1000   # scale for the energy colomap 
emax =  1000

#----------------

read_psi = LegacyVTKReader()
read_ener = LegacyVTKReader()

read_psi.FileNames = file_psi
read_ener.FileNames = file_ener

reader = AppendAttributes(read_psi, read_ener)
surface = ExtractSurface(reader)

psiwarp = WarpByScalar(surface)
psiwarp.Scalars="psi"
psiwarp.ScaleFactor=hscale

psiwarp_rep = Show(psiwarp)

lt = MakeBlueToRedLT(emin, emax)
lt.ColorSpace='Diverging'
psiwarp_rep.ColorArrayName="ener"
psiwarp_rep.ColorAttributeType=0    # 'Point_Data'?
psiwarp_rep.LookupTable=lt


Render()
view = GetActiveView()


view.CameraPosition=[2000, 511.5, 1000]
view.CameraViewUp = [-0.7, 0, 0.7]
view.CameraFocalPoint=[511.5, 511.5, 0]
view.CenterOfRotation=[511.5, 511.5, 0]

view.LightSwitch=0
view.UseLight=1

view.Background=[0.35, 0.35, 0.4]
view.KeyLightWarmth=0.6
view.KeyLightIntensity=0.75
view.KeyLightElevation=50
view.KeyLightAzimuth=10

view.FillLightWarmth=0.4
view.FillLightElevation=-75
view.FillLightAzimuth=-10
view.SetPropertyWithName("FillLightK:FRatio", 3.0)

view.BackLightWarmth=0.5
view.BackLightElevation=0
view.BackLightAzimuth=110
view.SetPropertyWithName("BackLightK:BRatio", 3.5)

view.HeadLightWarmth=0.5
view.SetPropertyWithName("HeadLightK:HRatio", 3.0)

view.ViewSize=(1200,800)

#-----------------------------------------


s =  "%04d" % n

file_img =  "tmp/surf" + s + ".jpg"

view.ViewTime=n
Render()
WriteImage(file_img)


#-----------------------------------------

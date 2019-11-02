
from paraview.simple import *


def process_one(n):

#-- parameters --

    hscale = 10.0  # scale fot the height of psi

    emin = -1000   # scale for the energy colomap 
    emax =  1000
    
    s =  "%04d" % n

    file_psi  = "eps1e-2_4096_zoom_VTK/psi." + s + ".vtk"
    file_ener = "eps1e-2_4096_zoom_VTK/ener." + s + ".vtk"
    file_img =  "IMG_surf_ener/surf" + s + ".jpg"

#----------------

    read_psi = LegacyVTKReader(FileNames=file_psi)
    read_ener = LegacyVTKReader(FileNames=file_ener)
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

    view = Render()
    view.CameraPosition=[1600, 511.5, 1600]
    view.CameraViewUp = [-0.7, 0, 0.7]
    view.CameraFocalPoint=[511.5, 511.5, -100]
    view.CenterOfRotation=[511.5, 511.5, -100]

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

    view.ViewSize=(800,600)

    WriteImage(file_img)

#-----------------------------------------

n0 = 47
nmax = 1120

for n in range(n0, 2+nmax, 1):
    process_one(n)

#-----------------------------------------

from paraview.simple import *
from math import *

# run as a batch: /Applications/paraview.app/Contents/bin/pvbatch guide_nkflow.py

#-- file selection --

#file_dir = "/Users/nata/projects/tmp/paraview/t15a20_"
#file_dir = "/Users/nata/projects/tmp/paraview/t25a14_"
file_dir = "/Users/nata/projects/tmp/paraview/t75a12_"
file_pk  = file_dir + "pk.vtk"
file_qk  = file_dir + "qk.vtk"
file_img = file_dir + "pk.png"


qmax = 0.001 # max charge, for colormap scale

#-- read potential and compute stream lines ---

charge    = LegacyVTKReader( FileNames=file_qk )
potential = LegacyVTKReader( FileNames=file_pk )

grad      = Gradient()
grad.SelectInputScalars = ['POINTS', 'pk']

r=4

#nmax=18;  a0=0.30               #t15a20           
#nmax=18;  a0=0.10 + pi/18/2     #t25a14
nmax=20;  a0=0.15 + pi/12       #t75a12 

for n in range(0, nmax, 1):
    
    a = a0 + 2*pi*n/nmax
    x = 256 + r*cos(a)
    y = 256 + r*sin(a)

    stream = StreamTracer(grad, SeedType="Point Source" )
    stream.SeedType.Radius = 0
    stream.SeedType.NumberOfPoints = 1
    stream.IntegrationDirection = 'FORWARD'
    stream.MaximumStreamlineLength = 200.0
    stream.Vectors = ['POINTS', 'pkGradient']
    stream.SeedType.Center = [x, y, 0.0]

    stream.SMProxy.InvokeEvent('UserEvent', 'HideWidget')

    streamRep = Show()
    streamRep.DiffuseColor = [0.0, 0.0, 0.0]
    #streamRep.DiffuseColor = [0.0, 0.0, 0.8]
    streamRep.ColorArrayName = ''
    streamRep.LineWidth=2

#-- charge representations --

#myRGB = [-2.0*qmax, 0.0, 0.0, 0.0, qmax, 1.0, 1.0, 1.0]
#lt = CreateLookupTable(RGBPoints=myRGB)

lt = MakeBlueToRedLT(-qmax, qmax)
lt.ColorSpace='Diverging'
lt.LockScalarRange = 1


SetActiveSource(charge)

chargeRep = Show()
chargeRep.ColorArrayName = 'qk'
chargeRep.ColorAttributeType = 'POINT_DATA'
chargeRep.LookupTable = lt


#-- create camera view --

Render()

view = GetRenderView()
view.OrientationAxesVisibility = 0
view.CenterAxesVisibility = 0

view.CameraParallelScale = 180
view.CameraPosition = [256.5, 256.5, 1024]
view.CameraFocalPoint = [256.5, 256.5, 0.0]
view.CenterOfRotation = [256.5, 256.5, 0.0]
view.CameraClippingRange = [512, 512]
view.ViewSize=(512,512)

view.CameraClippingRange = [512, 512]

Render()
#WriteImage(file_img, Magnification=2)



#-----------------------------------------

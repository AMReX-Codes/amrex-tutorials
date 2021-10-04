# script-version: 2.0
# Catalyst state generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [800, 800]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5000000132713467, 0.5000000132713467, 0.4984375119674951]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-1.1558770108238157, -1.7832016892543263, 2.568577146937055]
renderView1.CameraFocalPoint = [0.3924535407483471, 0.3045464097802382, 0.5444715088099857]
renderView1.CameraViewUp = [0.3582901498357902, 0.49918337862420037, 0.7889512805211583]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8526515725882607

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(800, 800)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML UniformGrid AMR Reader'
mesh_000006vth = XMLUniformGridAMRReader(registrationName='mesh')
mesh_000006vth.TimeArray = 'None'

# create a new 'XML PolyData Reader'
particles_000000_000006vtp = XMLPolyDataReader(registrationName='particles')
particles_000000_000006vtp.PointArrayStatus = ['u']
particles_000000_000006vtp.TimeArray = 'None'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from particles_000000_000006vtp
particles_000000_000006vtpDisplay = Show(particles_000000_000006vtp, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.RGBPoints = [1.000000007718748, 0.231373, 0.298039, 0.752941, 1.2045469413608973, 0.865003, 0.865003, 0.865003, 1.4090938750030466, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
particles_000000_000006vtpDisplay.Representation = 'Surface'
particles_000000_000006vtpDisplay.ColorArrayName = ['POINTS', 'u']
particles_000000_000006vtpDisplay.LookupTable = uLUT
particles_000000_000006vtpDisplay.SelectTCoordArray = 'None'
particles_000000_000006vtpDisplay.SelectNormalArray = 'None'
particles_000000_000006vtpDisplay.SelectTangentArray = 'None'
particles_000000_000006vtpDisplay.OSPRayScaleArray = 'u'
particles_000000_000006vtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
particles_000000_000006vtpDisplay.SelectOrientationVectors = 'None'
particles_000000_000006vtpDisplay.ScaleFactor = 0.09846483203582466
particles_000000_000006vtpDisplay.SelectScaleArray = 'None'
particles_000000_000006vtpDisplay.GlyphType = 'Arrow'
particles_000000_000006vtpDisplay.GlyphTableIndexArray = 'None'
particles_000000_000006vtpDisplay.GaussianRadius = 0.004923241601791233
particles_000000_000006vtpDisplay.SetScaleArray = ['POINTS', 'u']
particles_000000_000006vtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
particles_000000_000006vtpDisplay.OpacityArray = ['POINTS', 'u']
particles_000000_000006vtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
particles_000000_000006vtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
particles_000000_000006vtpDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
particles_000000_000006vtpDisplay.ScaleTransferFunction.Points = [-0.9902609286592124, 0.0, 0.5, 0.0, 0.9926869306161271, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
particles_000000_000006vtpDisplay.OpacityTransferFunction.Points = [-0.9902609286592124, 0.0, 0.5, 0.0, 0.9926869306161271, 1.0, 0.5, 0.0]

# show data from mesh_000006vth
mesh_000006vthDisplay = Show(mesh_000006vth, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
mesh_000006vthDisplay.Representation = 'Outline'
mesh_000006vthDisplay.ColorArrayName = [None, '']
mesh_000006vthDisplay.LineWidth = 2.0
mesh_000006vthDisplay.SelectTCoordArray = 'None'
mesh_000006vthDisplay.SelectNormalArray = 'None'
mesh_000006vthDisplay.SelectTangentArray = 'None'
mesh_000006vthDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
mesh_000006vthDisplay.SelectOrientationVectors = 'None'
mesh_000006vthDisplay.ScaleFactor = 0.1
mesh_000006vthDisplay.SelectScaleArray = 'None'
mesh_000006vthDisplay.GlyphType = 'Arrow'
mesh_000006vthDisplay.GlyphTableIndexArray = 'None'
mesh_000006vthDisplay.GaussianRadius = 0.005
mesh_000006vthDisplay.SetScaleArray = [None, '']
mesh_000006vthDisplay.ScaleTransferFunction = 'PiecewiseFunction'
mesh_000006vthDisplay.OpacityArray = [None, '']
mesh_000006vthDisplay.OpacityTransferFunction = 'PiecewiseFunction'
mesh_000006vthDisplay.DataAxesGrid = 'GridAxesRepresentation'
mesh_000006vthDisplay.PolarAxes = 'PolarAxesRepresentation'
mesh_000006vthDisplay.ScalarOpacityUnitDistance = 0.027063293868263706

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Title = 'u'
uLUTColorBar.ComponentTitle = 'Magnitude'

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
particles_000000_000006vtpDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [1.000000007718748, 0.0, 0.5, 0.0, 1.4090938750030466, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'pv_render_particles_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [800, 800]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = './'
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)

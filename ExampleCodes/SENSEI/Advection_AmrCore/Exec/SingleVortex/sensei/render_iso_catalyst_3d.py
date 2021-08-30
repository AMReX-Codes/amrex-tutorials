# script-version: 2.0
# Catalyst state generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [800, 800]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [3.6901377539352986, 0.5341442836834782, 0.3594256077982131]
renderView1.CameraFocalPoint = [0.34407253898406637, 0.5341442836834782, 0.3594256077982131]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

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
mesh_000004vth = XMLUniformGridAMRReader(registrationName='mesh')
mesh_000004vth.DefaultNumberOfLevels = 10

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=mesh_000004vth)
cellDatatoPointData1.CellDataArraytoprocess = ['GhostType', 'phi', 'vtkGhostType']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=cellDatatoPointData1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5, 0.5, 0.5]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from cellDatatoPointData1
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'
cellDatatoPointData1Display.ColorArrayName = [None, '']
cellDatatoPointData1Display.LineWidth = 2.0
cellDatatoPointData1Display.OSPRayScaleArray = 'GhostType'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'GhostType'
cellDatatoPointData1Display.ScaleFactor = 0.1
cellDatatoPointData1Display.SelectScaleArray = 'GhostType'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'GhostType'
cellDatatoPointData1Display.GaussianRadius = 0.005
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'GhostType']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'GhostType']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.016960924803825335

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 8.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 8.0, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'phi'
phiLUT = GetColorTransferFunction('phi')
phiLUT.RGBPoints = [0.9998837073660926, 0.231373, 0.298039, 0.752941, 1.3647609434384869, 0.865003, 0.865003, 0.865003, 1.7296381795108813, 0.705882, 0.0156863, 0.14902]
phiLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Wireframe'
slice1Display.ColorArrayName = ['POINTS', 'phi']
slice1Display.LookupTable = phiLUT
slice1Display.OSPRayScaleArray = 'GhostType'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'GhostType'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'GhostType'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'GhostType'
slice1Display.GaussianRadius = 0.005
slice1Display.SetScaleArray = ['POINTS', 'GhostType']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'GhostType']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 6.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 6.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for phiLUT in view renderView1
phiLUTColorBar = GetScalarBar(phiLUT, renderView1)
phiLUTColorBar.Orientation = 'Horizontal'
phiLUTColorBar.WindowLocation = 'AnyLocation'
phiLUTColorBar.Position = [0.4079687499999999, 0.09624999999999999]
phiLUTColorBar.Title = 'phi'
phiLUTColorBar.ComponentTitle = ''
phiLUTColorBar.ScalarBarLength = 0.33000000000000007

# set color bar visibility
phiLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'phi'
phiPWF = GetOpacityTransferFunction('phi')
phiPWF.Points = [0.9998837073660926, 0.0, 0.5, 0.0, 1.7296381795108813, 1.0, 0.5, 0.0]
phiPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'pv_image_3d_%.6ts.png'
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

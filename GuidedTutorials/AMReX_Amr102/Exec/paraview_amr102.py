# based on traces generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

import subprocess
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--frame_rate', type=int, default=15, help="Frame rate for generating movies, i.e. number of plots per second in the movie.")
parser.add_argument('-r', '--resolution', type=int, default=1024, help="(Square) resolution of output movie.")
parser.add_argument('-d', '--spacedim', type=int, default=3, help="Dimensionality of the problem: 2 or 3")
parser.add_argument('-w', '--shade_weights', type=int, default=0, help="Toggle on (1) or off (0) shading particles by their weights.")
args = parser.parse_args()

def generate_movie_3D(AllPlotFiles):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
    plt00.CellArrayStatus = []

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Properties modified on plt00
    plt00.CellArrayStatus = ['phi', 'proc', 'vfrac', 'xvel', 'yvel', 'zvel']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1166, 1176]
    renderView1.ViewSize = [1200, 1200]

    # get layout
    layout1 = GetLayout()

    # show data in view
    plt00Display = Show(plt00, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00Display.Representation = 'Outline'
    plt00Display.ColorArrayName = [None, '']
    plt00Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00Display.SelectOrientationVectors = 'None'
    plt00Display.ScaleFactor = 0.1
    plt00Display.SelectScaleArray = 'None'
    plt00Display.GlyphType = 'Arrow'
    plt00Display.GlyphTableIndexArray = 'None'
    plt00Display.GaussianRadius = 0.005
    plt00Display.SetScaleArray = [None, '']
    plt00Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00Display.OpacityArray = [None, '']
    plt00Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00Display.PolarAxes = 'PolarAxesRepresentation'
    plt00Display.ScalarOpacityUnitDistance = 0.04436647145156464

    # reset view to fit data
    renderView1.ResetCamera()

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Slice'
    slice1 = Slice(Input=plt00)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.5, 0.5, 0.0625]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.0625]

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1.SliceType)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.1
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 0.005
    slice1Display.SetScaleArray = [None, '']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = [None, '']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(slice1Display, ('FIELD', 'vtkBlockColors'))

    # get color transfer function/color map for 'vtkBlockColors'
    vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

    # get opacity transfer function/opacity map for 'vtkBlockColors'
    vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

    # set scalar coloring
    ColorBy(slice1Display, ('CELLS', 'phi'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'phi'
    phiLUT = GetColorTransferFunction('phi')

    # get opacity transfer function/opacity map for 'phi'
    phiPWF = GetOpacityTransferFunction('phi')

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get color legend/bar for phiLUT in view renderView1
    phiLUTColorBar = GetScalarBar(phiLUT, renderView1)

    # change scalar bar placement
    phiLUTColorBar.WindowLocation = 'AnyLocation'
    phiLUTColorBar.Position = [0, 0.75]
    phiLUTColorBar.ScalarBarLength = 0.2

    # create a new 'AMReX/BoxLib Particles Reader'
    plt00_1 = AMReXBoxLibParticlesReader(FileNames=AllPlotFiles)
    plt00_1.PointArrayStatus = ['id', 'cpu', 'real_comp0', 'real_comp1', 'real_comp2', 'real_comp3']

    # show data in view
    plt00_1Display = Show(plt00_1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    plt00_1Display.Representation = 'Surface'
    plt00_1Display.ColorArrayName = [None, '']
    plt00_1Display.OSPRayScaleArray = 'cpu'
    plt00_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00_1Display.SelectOrientationVectors = 'None'
    plt00_1Display.ScaleFactor = 0.029962979569768036
    plt00_1Display.SelectScaleArray = 'None'
    plt00_1Display.GlyphType = 'Arrow'
    plt00_1Display.GlyphTableIndexArray = 'None'
    plt00_1Display.GaussianRadius = 0.0014981489784884018
    plt00_1Display.SetScaleArray = ['POINTS', 'cpu']
    plt00_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00_1Display.OpacityArray = ['POINTS', 'cpu']
    plt00_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00_1Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00_1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    plt00_1Display.ScaleTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    plt00_1Display.OpacityTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(plt00_1Display, ('FIELD', 'vtkBlockColors'))

    # create a new 'Glyph'
    glyph1 = Glyph(Input=plt00_1,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'No orientation array']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 0.029962979569768036
    glyph1.GlyphTransform = 'Transform2'

    # Properties modified on glyph1
    glyph1.GlyphType = 'Sphere'
    glyph1.ScaleFactor = 0.01
    glyph1.GlyphMode = 'Every Nth Point'
    glyph1.Stride = 100

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'Normals'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.030791985988616946
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 0.0015395992994308473
    glyph1Display.SetScaleArray = ['POINTS', 'Normals']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'Normals']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    glyph1Display.ScaleTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    glyph1Display.OpacityTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(glyph1Display, ('FIELD', 'vtkBlockColors'))

    if args.shade_weights:
        # set scalar coloring
        ColorBy(glyph1Display, ('POINTS', 'real_comp3'))

        # rescale color and/or opacity maps used to include current data range
        glyph1Display.RescaleTransferFunctionToDataRange(True, False)

        # get color transfer function/color map for 'real_comp3'
        real_comp3LUT = GetColorTransferFunction('real_comp3')

        # get opacity transfer function/opacity map for 'real_comp3'
        real_comp3PWF = GetOpacityTransferFunction('real_comp3')

        # show color bar/color legend
        glyph1Display.SetScalarBarVisibility(renderView1, True)

        # get color legend/bar for real_comp3LUT in view renderView1
        real_comp3LUTColorBar = GetScalarBar(real_comp3LUT, renderView1)

        # change scalar bar placement
        real_comp3LUTColorBar.WindowLocation = 'AnyLocation'
        real_comp3LUTColorBar.Position = [0, 0.4]
        real_comp3LUTColorBar.ScalarBarLength = 0.2

    # create a new 'XML Partitioned Polydata Reader'
    ebpvtp = XMLPartitionedPolydataReader(FileName=['eb.pvtp'])

    # show data in view
    ebpvtpDisplay = Show(ebpvtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    ebpvtpDisplay.Representation = 'Surface'
    ebpvtpDisplay.ColorArrayName = [None, '']
    ebpvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    ebpvtpDisplay.SelectOrientationVectors = 'None'
    ebpvtpDisplay.ScaleFactor = 0.019999998807907107
    ebpvtpDisplay.SelectScaleArray = 'None'
    ebpvtpDisplay.GlyphType = 'Arrow'
    ebpvtpDisplay.GlyphTableIndexArray = 'None'
    ebpvtpDisplay.GaussianRadius = 0.0009999999403953552
    ebpvtpDisplay.SetScaleArray = [None, '']
    ebpvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.OpacityArray = [None, '']
    ebpvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
    ebpvtpDisplay.PolarAxes = 'PolarAxesRepresentation'

    # update the view to ensure updated data information
    renderView1.Update()

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.806487938976244, 0.6098478265765159, 2.3240231832552687]
    renderView1.CameraFocalPoint = [0.4999999999999996, 0.49999999999999956, 0.06249999999999997]
    renderView1.CameraViewUp = [-0.004859468445493335, 0.9988423431014976, -0.04785769733216954]
    renderView1.CameraParallelScale = 0.7159515667518355

    # save animation
    output_movie_base = "amr102_3D"
    output_movie = output_movie_base + ".avi"
    SaveAnimation(output_movie,
                  renderView1,
                  ImageResolution=[1200, 1200],
                  FrameRate=args.frame_rate,
                  FrameWindow=[0, len(AllPlotFiles)-1])

    return output_movie_base, output_movie

def generate_movie_2D(AllPlotFiles):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
    plt00.CellArrayStatus = []

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Properties modified on plt00
    plt00.CellArrayStatus = ['phi', 'proc', 'vfrac', 'xvel', 'yvel']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1166, 1176]
    renderView1.ViewSize = [1200, 1200]

    # get layout
    layout1 = GetLayout()

    # show data in view
    plt00Display = Show(plt00, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00Display.Representation = 'Outline'
    plt00Display.ColorArrayName = [None, '']
    plt00Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00Display.SelectOrientationVectors = 'None'
    plt00Display.ScaleFactor = 0.1
    plt00Display.SelectScaleArray = 'None'
    plt00Display.GlyphType = 'Arrow'
    plt00Display.GlyphTableIndexArray = 'None'
    plt00Display.GaussianRadius = 0.005
    plt00Display.SetScaleArray = [None, '']
    plt00Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00Display.OpacityArray = [None, '']
    plt00Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00Display.PolarAxes = 'PolarAxesRepresentation'
    plt00Display.ScalarOpacityUnitDistance = 0.08838834764831845

    # reset view to fit data
    renderView1.ResetCamera()

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # change representation type
    plt00Display.SetRepresentationType('Surface')

    # set scalar coloring
    ColorBy(plt00Display, ('CELLS', 'phi'))

    # rescale color and/or opacity maps used to include current data range
    plt00Display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'phi'
    phiLUT = GetColorTransferFunction('phi')

    # get opacity transfer function/opacity map for 'phi'
    phiPWF = GetOpacityTransferFunction('phi')

    # show color bar/color legend
    plt00Display.SetScalarBarVisibility(renderView1, True)

    # get color legend/bar for phiLUT in view renderView1
    phiLUTColorBar = GetScalarBar(phiLUT, renderView1)

    # change scalar bar placement
    phiLUTColorBar.WindowLocation = 'AnyLocation'
    phiLUTColorBar.Position = [0, 0.75]
    phiLUTColorBar.ScalarBarLength = 0.2

    # create a new 'AMReX/BoxLib Particles Reader'
    plt00_1 = AMReXBoxLibParticlesReader(FileNames=AllPlotFiles)
    plt00_1.PointArrayStatus = ['id', 'cpu', 'real_comp0', 'real_comp1', 'real_comp2']

    # show data in view
    plt00_1Display = Show(plt00_1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    plt00_1Display.Representation = 'Surface'
    plt00_1Display.ColorArrayName = [None, '']
    plt00_1Display.OSPRayScaleArray = 'cpu'
    plt00_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00_1Display.SelectOrientationVectors = 'None'
    plt00_1Display.ScaleFactor = 0.08994852971110677
    plt00_1Display.SelectScaleArray = 'None'
    plt00_1Display.GlyphType = 'Arrow'
    plt00_1Display.GlyphTableIndexArray = 'None'
    plt00_1Display.GaussianRadius = 0.004497426485555338
    plt00_1Display.SetScaleArray = ['POINTS', 'cpu']
    plt00_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00_1Display.OpacityArray = ['POINTS', 'cpu']
    plt00_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00_1Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00_1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    plt00_1Display.ScaleTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    plt00_1Display.OpacityTransferFunction.Points = [2.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(plt00_1Display, ('FIELD', 'vtkBlockColors'))

    # get color transfer function/color map for 'vtkBlockColors'
    vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

    # get opacity transfer function/opacity map for 'vtkBlockColors'
    vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

    # create a new 'Glyph'
    glyph1 = Glyph(Input=plt00_1,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'No orientation array']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 0.08994852971110677
    glyph1.GlyphTransform = 'Transform2'

    # Properties modified on glyph1
    glyph1.GlyphType = 'Sphere'
    glyph1.ScaleFactor = 0.01
    glyph1.GlyphMode = 'Every Nth Point'
    glyph1.Stride = 100

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'Normals'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.0902662232572228
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 0.00451331116286114
    glyph1Display.SetScaleArray = ['POINTS', 'Normals']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'Normals']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    glyph1Display.ScaleTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    glyph1Display.OpacityTransferFunction.Points = [-0.9749279618263245, 0.0, 0.5, 0.0, 0.9749279618263245, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(glyph1Display, ('FIELD', 'vtkBlockColors'))

    if args.shade_weights:
        # set scalar coloring
        ColorBy(glyph1Display, ('POINTS', 'real_comp2'))

        # rescale color and/or opacity maps used to include current data range
        glyph1Display.RescaleTransferFunctionToDataRange(True, False)

        # get color transfer function/color map for 'real_comp2'
        real_comp2LUT = GetColorTransferFunction('real_comp2')

        # get opacity transfer function/opacity map for 'real_comp2'
        real_comp2PWF = GetOpacityTransferFunction('real_comp2')

        # show color bar/color legend
        glyph1Display.SetScalarBarVisibility(renderView1, True)

        # get color legend/bar for real_comp2LUT in view renderView1
        real_comp2LUTColorBar = GetScalarBar(real_comp2LUT, renderView1)

        # change scalar bar placement
        real_comp2LUTColorBar.WindowLocation = 'AnyLocation'
        real_comp2LUTColorBar.Position = [0, 0.4]
        real_comp2LUTColorBar.ScalarBarLength = 0.2

    # create a new 'XML Partitioned Polydata Reader'
    ebpvtp = XMLPartitionedPolydataReader(FileName=['eb.pvtp'])

    # show data in view
    ebpvtpDisplay = Show(ebpvtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    ebpvtpDisplay.Representation = 'Surface'
    ebpvtpDisplay.ColorArrayName = [None, '']
    ebpvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    ebpvtpDisplay.SelectOrientationVectors = 'None'
    ebpvtpDisplay.ScaleFactor = 0.019999998807907107
    ebpvtpDisplay.SelectScaleArray = 'None'
    ebpvtpDisplay.GlyphType = 'Arrow'
    ebpvtpDisplay.GlyphTableIndexArray = 'None'
    ebpvtpDisplay.GaussianRadius = 0.0009999999403953552
    ebpvtpDisplay.SetScaleArray = [None, '']
    ebpvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.OpacityArray = [None, '']
    ebpvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
    ebpvtpDisplay.PolarAxes = 'PolarAxesRepresentation'

    # update the view to ensure updated data information
    renderView1.Update()

    # change representation type
    ebpvtpDisplay.SetRepresentationType('Wireframe')

    # change solid color
    ebpvtpDisplay.AmbientColor = [0.0, 0.0, 0.0]
    ebpvtpDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # Properties modified on ebpvtpDisplay
    ebpvtpDisplay.LineWidth = 2.0

    # Properties modified on ebpvtpDisplay
    ebpvtpDisplay.LineWidth = 3.0

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = 0.5893976543919168

    # save animation
    output_movie_base = "amr102_2D"
    output_movie = output_movie_base + ".avi"
    SaveAnimation(output_movie,
                  renderView1,
                  ImageResolution=[1200, 1200],
                  FrameRate=args.frame_rate,
                  FrameWindow=[0, len(AllPlotFiles)-1])

    return output_movie_base, output_movie

def convert_avi_to_gif(output_movie_base, output_movie):
    # use ffmpeg to convert the avi movie into an animated gif
    ffmpeg_convert_to_gif = 'ffmpeg -y -i {} -vf "fps=35,scale={}:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 {}.gif'.format(output_movie, args.resolution, output_movie_base)
    subprocess.run(ffmpeg_convert_to_gif, shell=True)

if __name__ == "__main__":
    if not (args.spacedim == 2 or args.spacedim == 3):
        print("Please specify --spacedim D (with D=2 or D=3)")
        exit()

    if args.frame_rate <= 0:
        print("Please specify --frame_rate F (with F > 0)")
        exit()

    if args.resolution <= 0:
        print("Please specify --resolution R (with R > 0)")
        exit()

    if not (args.shade_weights == 0 or args.shade_weights == 1):
        print("Please specify --shade_weights S (with S = 0 or S = 1)")
        exit()

    # get all the plotfiles
    PlotFiles = sorted(glob.glob("plt" + "[0-9]"*5))

    # call the 2D or 3D vis function
    output_movie_base = None
    output_movie = None

    if args.spacedim == 3:
        output_movie_base, output_movie = generate_movie_3D(PlotFiles)
    elif args.spacedim == 2:
        output_movie_base, output_movie = generate_movie_2D(PlotFiles)

    # convert the avi movie into an animated gif
    convert_avi_to_gif(output_movie_base, output_movie)
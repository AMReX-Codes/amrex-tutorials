# trace generated using paraview version 5.8.0
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
args = parser.parse_args()

def generate_movie_3D(AllPlotFiles):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'XML Partitioned Polydata Reader'
    ebpvtp = XMLPartitionedPolydataReader(FileName=['eb.pvtp'])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1200, 1200]

    # get layout
    layout1 = GetLayout()

    # show data in view
    ebpvtpDisplay = Show(ebpvtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    ebpvtpDisplay.Representation = 'Surface'
    ebpvtpDisplay.ColorArrayName = [None, '']
    ebpvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    ebpvtpDisplay.SelectOrientationVectors = 'None'
    ebpvtpDisplay.ScaleFactor = 0.14000000208616256
    ebpvtpDisplay.SelectScaleArray = 'None'
    ebpvtpDisplay.GlyphType = 'Arrow'
    ebpvtpDisplay.GlyphTableIndexArray = 'None'
    ebpvtpDisplay.GaussianRadius = 0.007000000104308128
    ebpvtpDisplay.SetScaleArray = [None, '']
    ebpvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.OpacityArray = [None, '']
    ebpvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    ebpvtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
    ebpvtpDisplay.PolarAxes = 'PolarAxesRepresentation'

    # reset view to fit data
    renderView1.ResetCamera()

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00000 = AMReXBoxLibGridReader(FileNames=['plt00000'])
    plt00000.CellArrayStatus = []

    # Properties modified on plt00000
    plt00000.CellArrayStatus = ['proc', 'vfrac']

    # show data in view
    plt00000Display = Show(plt00000, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00000Display.Representation = 'Outline'
    plt00000Display.ColorArrayName = [None, '']
    plt00000Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00000Display.SelectOrientationVectors = 'None'
    plt00000Display.ScaleFactor = 0.2
    plt00000Display.SelectScaleArray = 'None'
    plt00000Display.GlyphType = 'Arrow'
    plt00000Display.GlyphTableIndexArray = 'None'
    plt00000Display.GaussianRadius = 0.01
    plt00000Display.SetScaleArray = [None, '']
    plt00000Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00000Display.OpacityArray = [None, '']
    plt00000Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00000Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00000Display.PolarAxes = 'PolarAxesRepresentation'
    plt00000Display.ScalarOpacityUnitDistance = 0.051941539345088994

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Slice'
    slice1 = Slice(Input=plt00000)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.625, 1.0, 0.25]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.625, 1.0, 0.25]

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
    slice1Display.ScaleFactor = 0.2
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 0.01
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
    ColorBy(slice1Display, ('CELLS', 'proc'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'proc'
    procLUT = GetColorTransferFunction('proc')

    # get opacity transfer function/opacity map for 'proc'
    procPWF = GetOpacityTransferFunction('proc')

    # create a new 'AMReX/BoxLib Particles Reader'
    plt0 = AMReXBoxLibParticlesReader(FileNames=AllPlotFiles)
    plt0.PointArrayStatus = ['id', 'cpu', 'int_comp0', 'real_comp0', 'real_comp1', 'real_comp2']

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # show data in view
    plt0Display = Show(plt0, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    plt0Display.Representation = 'Surface'
    plt0Display.ColorArrayName = [None, '']
    plt0Display.OSPRayScaleArray = 'cpu'
    plt0Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt0Display.SelectOrientationVectors = 'None'
    plt0Display.ScaleFactor = 0.11000000000000001
    plt0Display.SelectScaleArray = 'None'
    plt0Display.GlyphType = 'Arrow'
    plt0Display.GlyphTableIndexArray = 'None'
    plt0Display.GaussianRadius = 0.0055000000000000005
    plt0Display.SetScaleArray = ['POINTS', 'cpu']
    plt0Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt0Display.OpacityArray = ['POINTS', 'cpu']
    plt0Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt0Display.DataAxesGrid = 'GridAxesRepresentation'
    plt0Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    plt0Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    plt0Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(plt0Display, ('FIELD', 'vtkBlockColors'))

    # create a new 'Glyph'
    glyph1 = Glyph(Input=plt0,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'No orientation array']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 0.11000000000000001
    glyph1.GlyphTransform = 'Transform2'

    # Properties modified on glyph1
    glyph1.GlyphType = 'Sphere'
    glyph1.ScaleFactor = 0.05
    glyph1.GlyphMode = 'All Points'

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'Normals'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.11487463433295489
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 0.005743731716647744
    glyph1Display.SetScaleArray = ['POINTS', 'Normals']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'Normals']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    glyph1Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    glyph1Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(glyph1Display, ('FIELD', 'vtkBlockColors'))

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5999999884516001, 0.9000000134110451, 5.480672965279949]
    renderView1.CameraFocalPoint = [0.5999999884516001, 0.9000000134110451, 0.25]
    renderView1.CameraParallelScale = 0.9246621010295244

    # save animation
    output_movie_base = "pachinko"
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
    if args.frame_rate <= 0:
        print("Please specify --frame_rate F (with F > 0)")
        exit()

    if args.resolution <= 0:
        print("Please specify --resolution R (with R > 0)")
        exit()

    # get all the plotfiles
    PlotFiles = sorted(glob.glob("plt" + "[0-9]"*5))

    # generate the avi movie file
    output_movie_base, output_movie = generate_movie_3D(PlotFiles)

    # convert the avi movie into an animated gif
    convert_avi_to_gif(output_movie_base, output_movie)
# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2019 replay file
# Internal Version: 2018_09_24-13.41.51 157541
# Run by nidish on Wed Mar 18 10:18:44 2020
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=401.108306884766, 
    height=160.477783203125)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.ModelFromInputFile(name='singbrick_cons', 
    inputFileName='/home/nidish/Documents/Academics/RESEARCH/INTRED_STUDY/ABAQUS_MODELS/SINGBRICK/singbrick_cons.inp')
#: The model "singbrick_cons" has been created.
#: The part "PART-1" has been imported from the input file.
#: 
#: WARNING: The following keywords/parameters are not yet supported by the input file reader:
#: ---------------------------------------------------------------------------------
#: *MATRIXGENERATE
#: *MATRIXOUTPUT
#: *PREPRINT
#: The model "singbrick_cons" has been imported from an input file. 
#: Please scroll up to check for error and warning messages.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['singbrick_cons'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
    interactions=ON, constraints=ON, connectors=ON, engineeringFeatures=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.50246, 
    farPlane=9.25125, width=5.4242, height=2.5821, cameraPosition=(4.21373, 
    -0.193559, 6.01932), cameraUpVector=(-0.189261, 0.947882, -0.256321), 
    cameraTarget=(-0.00320324, -0.0791156, 0.0406874))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.55001, 
    farPlane=9.21801, width=5.47107, height=2.60441, cameraPosition=(1.96449, 
    -1.11555, 7.01561), cameraUpVector=(0.0683435, 0.971714, -0.226055), 
    cameraTarget=(-0.0214298, -0.0865869, 0.0487608))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.22889, 
    farPlane=9.6195, width=5.15452, height=2.45372, cameraPosition=(4.86472, 
    -4.01131, 3.85722), cameraUpVector=(-0.195115, 0.855132, 0.480291), 
    cameraTarget=(0.00485795, -0.112834, 0.0201331))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.46701, 
    farPlane=9.3994, width=5.38926, height=2.56547, cameraPosition=(6.13308, 
    -3.97812, 1.10482), cameraUpVector=(-0.0478985, 0.747824, 0.662167), 
    cameraTarget=(0.0231585, -0.112355, -0.0195796))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.68407, 
    farPlane=9.10945, width=5.60324, height=2.66733, cameraPosition=(7.08569, 
    -0.165543, 1.9578), cameraUpVector=(-0.459299, 0.765722, 0.450237), 
    cameraTarget=(0.0380413, -0.0527904, -0.00625328))
session.viewports['Viewport: 1'].view.setValues(nearPlane=5.38878, 
    farPlane=9.35675, width=5.31215, height=2.52876, cameraPosition=(4.23723, 
    0.790579, 5.94891), cameraUpVector=(-0.297012, 0.897324, -0.326485), 
    cameraTarget=(0.00735399, -0.0424898, 0.0367441))

# OpenVDBsimulations

Overview of what files contain:

Helper Files:

ppm.h - includes methods for displaying .ppm files (images)

STLParser.h - parses binary STL files for the OpenVDB format

dataConverter.cpp - extracts particle data from .vdb files into .csv files for easier post-processing

Basic Files:

PointCload.cpp - generates a bunch of random points

Points.cpp - example of getting particle properties

ParticleViewer.cpp - example of rendering particles as level-sets

Early Simulations:

  These are some initial simulations that were performed
  Some of the methods used in these simulations are not used in later simulations
  becuase they take too long but may be revisited later

Stationary.cpp - example of stationary sphere with particles moving around it

Moving.cpp - example of stationary particles with sphere moving throught them

Rotation.cpp - example of cube rotationg while moving through stationary particles

STLobject.cpp - example of a moving object, where the object is specified by the user at the command line with an STL file name

cube.stl is an stl of a cube

disks.stl is an stl of 2 disks at 90 degress to each other
  
Later Simulations:

diverseTest.cpp - example of a cough with different particle densities displayed in different colors (blue are smaller, red is bigger)

cubeTest.cpp - same as diverseTest, but with a cube moving through the coughed particles


# ImageSpaceVis

The repsoitory contains updated version of the code that can be compiled in Microsoft Visual Studio 2019.

## Pre-Requisites:

- OpenGL
- CMake GUI

## Build Steps:

- Download or clone the repository
- Make a new directory named "build" in the main directory of the project
- Open CMake GUI
- Set the source code path to project directory
- Set the build path to build directory previously made in the in the project directory
- Configure (Make sure to set the generator to Visual Studio 2019 and platform for generator to Win32)
- Generate
- Open the Project (.sln file in build directory) using Visual Studio
- Build the project

### Additional Step:

- Add the GL/lib path found in project directory to the environment PATH variable to make the shared libraries (.dll) available to the build binaries

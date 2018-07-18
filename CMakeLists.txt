cmake_minimum_required(VERSION 2.8.9)

project(NeuroPath)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/NeuroPath")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Seth Baksman (UWO)")
set(EXTENSION_DESCRIPTION "A 3D Slicer Widget for Neurosurgical Trajectory Visualization")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/NeuroPath.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/NeuroPath/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(NeuroPath)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})

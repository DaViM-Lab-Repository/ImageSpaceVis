cmake_minimum_required (VERSION 2.8)

project (ImageSpaceVis)
# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
  SET (CMAKE_BUILD_TYPE Release)
ENDIF()

add_subdirectory(src)





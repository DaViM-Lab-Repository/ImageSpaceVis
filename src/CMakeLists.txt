set(LIB_DIR ${PROJECT_SOURCE_DIR}/GL/lib)
include_directories(
    ${PROJECT_SOURCE_DIR}/GL/include
)
link_directories(
    ${LIB_DIR}
)

# find_library(OPENGL OpenGL32 PATHS ${PROJECT_SOURCE_DIR}/GL/lib)
find_library(GLUT glut32 PATHS ${PROJECT_SOURCE_DIR}/GL/lib)
find_library(GLUI glui32 PATHS ${PROJECT_SOURCE_DIR}/GL/lib)
# find_library(GLU glu32 PATHS ${PROJECT_SOURCE_DIR}/GL/lib)
# find_library(GLAUX glaux PATHS ${PROJECT_SOURCE_DIR}/GL/lib)
add_library(app SHARED STATIC
    directed_graph.h
    directed_graph.cpp
    Geometry.cpp
    icMatrix.h
    icVector.h
    icVector.cpp
    Img_MorseDecomp.h
    Img_MorseDecomp.cpp
    ply.h
    ply.cpp
    Skeleton.h
    VFAnalysis.h
    VFAnalysis.cpp
)
find_package(OpenGL REQUIRED)
if (OpenGL_FOUND)
    include_directories(${OpenGL_INCLUDE_DIRS})
    target_link_libraries(app ${OpenGL_LIBS})    
endif()
target_link_libraries(app ${GLUT} ${GLUI})
add_executable(sample sample.cpp)
target_link_libraries(sample app)
cmake_minimum_required (VERSION 2.6)
project (FMMMoM)

# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)

set (CMAKE_CXX_FLAGS "-g -Wall -std=c++11")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/FMMMoM_Config.h.in"
  "${PROJECT_BINARY_DIR}/include/FMMMoM_Config.h"
  )

# add the binary tree to the search path for include files
# so that we will find FMMMoM_Config.h
include_directories ("${PROJECT_BINARY_DIR}/include")

include_directories ("${PROJECT_SOURCE_DIR}/include")

file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# add the executable
add_executable (FieldData  ${SOURCES} ${PROJECT_SOURCE_DIR}/src/FieldData.cxx)
target_link_libraries(FieldData -larmadillo)
if(VTK_LIBRARIES)
  target_link_libraries(FieldData ${VTK_LIBRARIES})
else()
  target_link_libraries(FieldData vtkHybrid vtkWidgets)
endif()

# add the install targets
#install (TARGETS Tutorial DESTINATION bin)
#install (FILES "${PROJECT_BINARY_DIR}/include/FMMMoM_Config.h"
#  DESTINATION include)



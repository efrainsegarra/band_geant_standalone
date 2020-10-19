# Set project name
PROJECT(sigmaep_bound)
SET(LIBNAME "sigmaep_bound")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Cross Sections Library
add_library(${LIBNAME} sigmaep_bound.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES})

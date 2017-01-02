# CMakeLists for the library containing Wim's cross section codes

# Set up
cmake_minimum_required (VERSION 2.6)
enable_language(Fortran)

project (wim_lib)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

file( GLOB HEADERS ./*.h )

add_library(wim_lib crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f ${HEADERS} )

# install the grid files so that the cross section code can access them
INSTALL (FILES grids/a09.dsfs_lNNC grids/a09.sfs_lNNC
	DESTINATION "~/.deuteron_dis/grids/")

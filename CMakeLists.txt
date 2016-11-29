# CMakeLists for Deuteron DIS project

# Set up
cmake_minimum_required (VERSION 2.6)
enable_language(Fortran)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (deuteron_dis)

file( GLOB HEADERS ${CMAKE_SOURCE_DIR}/*.h )

set(INCLUDES ${ROOT_INCLUDE_DIR} )
include_directories(${INCLUDES})

add_executable(dis_CLAS12 dis_CLAS12.cpp crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f ${HEADERS} )
target_link_libraries(dis_CLAS12 ${ROOT_LIBRARIES} Foam)

add_executable(inclusive_CLAS12 inclusive_CLAS12.cpp crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f ${HEADERS} )

cmake_minimum_required (VERSION 2.6)
project (deuteron_dis)

enable_language(Fortran)

file( GLOB HEADERS ${CMAKE_SOURCE_DIR}/*.h )

add_executable(dis_CLAS12 dis_CLAS12.cpp crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f)

add_executable(inclusive_CLAS12 dis_CLAS12.cpp crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f)

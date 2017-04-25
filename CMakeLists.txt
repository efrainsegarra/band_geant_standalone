# CMakeLists for Deuteron DIS project

# Set up
cmake_minimum_required (VERSION 2.6)
enable_language(Fortran)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)
find_package(Geant4 REQUIRED ui_all vis_all)

project (deuteron_dis)

add_subdirectory(detector_sim)
add_subdirectory(wim_lib)
add_subdirectory(analysis)
add_subdirectory(tree)
add_subdirectory(lad_sim)

file( GLOB HEADERS ${CMAKE_SOURCE_DIR}/*.h )

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(dis_CLAS12 dis_CLAS12.cpp ${HEADERS} )
target_link_libraries(dis_CLAS12 ${ROOT_LIBRARIES} Foam wim_lib deuteron_dis_tree)

add_executable(dis_CLAS12_2d dis_CLAS12_2d.cpp ${HEADERS} )
target_link_libraries(dis_CLAS12_2d ${ROOT_LIBRARIES} Foam wim_lib deuteron_dis_tree)

add_executable(inclusive_CLAS12 inclusive_CLAS12.cpp ${HEADERS} )
target_link_libraries(inclusive_CLAS12 ${ROOT_LIBRARIES} Foam wim_lib deuteron_dis_tree)

add_executable(randoms_CLAS12 randoms_CLAS12.cpp )
target_link_libraries(randoms_CLAS12 ${ROOT_LIBRARIES} deuteron_dis_tree)

add_executable(isotropic_CLAS12 isotropic_CLAS12.cpp )
target_link_libraries(isotropic_CLAS12 ${ROOT_LIBRARIES} deuteron_dis_tree)

add_executable(dis_HallC dis_HallC.cpp ${HEADERS} )
target_link_libraries(dis_HallC ${ROOT_LIBRARIES} Foam wim_lib deuteron_dis_tree)

add_executable(inclusive_HallC inclusive_HallC.cpp ${HEADERS} )
target_link_libraries(inclusive_HallC ${ROOT_LIBRARIES} Foam wim_lib deuteron_dis_tree)

add_executable(randoms_HallC randoms_HallC.cpp )
target_link_libraries(randoms_HallC ${ROOT_LIBRARIES} deuteron_dis_tree)

cmake_minimum_required (VERSION 2.6)

project(lad_sim)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(lad_sim ${CMAKE_CURRENT_SOURCE_DIR}/lad_sim.cpp lad_helpers.cpp)
target_link_libraries(lad_sim ${ROOT_LIBRARIES} deuteron_dis_tree )

add_executable(lad_digi ${CMAKE_CURRENT_SOURCE_DIR}/lad_digi.cpp lad_helpers.cpp)
target_link_libraries(lad_digi ${ROOT_LIBRARIES} deuteron_dis_tree )

add_executable(lad_analysis ${CMAKE_CURRENT_SOURCE_DIR}/lad_analysis.cpp lad_helpers.cpp)
target_link_libraries(lad_analysis ${ROOT_LIBRARIES} deuteron_dis_tree )

add_executable(lad_optimize ${CMAKE_CURRENT_SOURCE_DIR}/lad_optimize.cpp lad_helpers.cpp)
target_link_libraries(lad_optimize ${ROOT_LIBRARIES} deuteron_dis_tree )



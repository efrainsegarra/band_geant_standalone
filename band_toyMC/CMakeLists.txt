cmake_minimum_required (VERSION 2.6)

project(band_toyMC)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(det_sim ${CMAKE_CURRENT_SOURCE_DIR}/det_sim.cpp)
target_link_libraries(det_sim ${ROOT_LIBRARIES} deuteron_dis_tree )

add_executable(make_res_hists ${CMAKE_CURRENT_SOURCE_DIR}/make_res_hists.cpp)
target_link_libraries(make_res_hists ${ROOT_LIBRARIES} )

cmake_minimum_required (VERSION 2.6)

#project(make_res_hists)

set(INCLUDES ${ROOT_INCLUDE_DIR} )
include_directories(${INCLUDES})

add_executable(make_res_hists ${CMAKE_CURRENT_SOURCE_DIR}/make_res_hists.cpp)
target_link_libraries(make_res_hists ${ROOT_LIBRARIES} )

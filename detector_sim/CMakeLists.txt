cmake_minimum_required (VERSION 2.6)

#project(det_sim)

set(INCLUDES ${ROOT_INCLUDE_DIR} )
include_directories(${INCLUDES})

add_executable(det_sim ${CMAKE_CURRENT_SOURCE_DIR}/det_sim.cpp)
target_link_libraries(det_sim ${ROOT_LIBRARIES} )

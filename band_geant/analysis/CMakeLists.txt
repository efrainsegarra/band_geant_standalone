cmake_minimum_required (VERSION 2.6)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
include_directories(${INCLUDES})

add_executable(makeHists ${CMAKE_CURRENT_SOURCE_DIR}/makeHists.cpp)
target_link_libraries(makeHists deuteron_dis_tree ${ROOT_LIBRARIES} )

add_executable(study_geantThreshold ${CMAKE_CURRENT_SOURCE_DIR}/study_geantThreshold.cpp )
target_link_libraries(study_geantThreshold deuteron_dis_tree ${ROOT_LIBRARIES})

add_executable(digitize ${CMAKE_CURRENT_SOURCE_DIR}/digitize.cpp )
target_link_libraries(digitize deuteron_dis_tree ${ROOT_LIBRARIES})

add_executable(resolutions ${CMAKE_CURRENT_SOURCE_DIR}/resolutions.cpp )
target_link_libraries(resolutions deuteron_dis_tree ${ROOT_LIBRARIES})
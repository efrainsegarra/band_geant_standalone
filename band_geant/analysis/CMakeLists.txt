cmake_minimum_required (VERSION 2.6)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
include_directories(${INCLUDES})

add_executable(peakToTail ${CMAKE_CURRENT_SOURCE_DIR}/peakToTail.cpp)
target_link_libraries(peakToTail deuteron_dis_tree ${ROOT_LIBRARIES} )

add_executable(study_geantThreshold ${CMAKE_CURRENT_SOURCE_DIR}/study_geantThreshold.cpp )
target_link_libraries(study_geantThreshold deuteron_dis_tree ${ROOT_LIBRARIES})

add_executable(digitize ${CMAKE_CURRENT_SOURCE_DIR}/digitize.cpp )
target_link_libraries(digitize deuteron_dis_tree ${ROOT_LIBRARIES})

add_executable(kinematicVars ${CMAKE_CURRENT_SOURCE_DIR}/kinematicVars.cpp )
target_link_libraries(kinematicVars deuteron_dis_tree ${ROOT_LIBRARIES})

add_executable(resHists ${CMAKE_CURRENT_SOURCE_DIR}/resHists.cpp)
target_link_libraries(resHists deuteron_dis_tree ${ROOT_LIBRARIES} )

add_executable(combineFiles ${CMAKE_CURRENT_SOURCE_DIR}/combineFiles.cpp)
target_link_libraries(combineFiles deuteron_dis_tree ${ROOT_LIBRARIES} )
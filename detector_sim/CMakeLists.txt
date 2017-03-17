cmake_minimum_required (VERSION 2.6)

project(det_sim)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories(${INCLUDES})

include(${Geant4_USE_FILE})

add_executable(det_sim ${CMAKE_CURRENT_SOURCE_DIR}/det_sim.cpp)
target_link_libraries(det_sim ${ROOT_LIBRARIES} deuteron_dis_tree )

file(GLOB G4sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB G4headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

add_executable(geant_sim ${CMAKE_CURRENT_SOURCE_DIR}/geant_sim.cpp ${G4sources} ${G4headers})
target_link_libraries(geant_sim ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} deuteron_dis_tree)




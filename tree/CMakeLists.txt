# Set project name
PROJECT(deuteron_dis_tree)

# find source and header files, set paths, includes...
SET(LIBNAME "deuteron_dis_tree")
FILE(GLOB SRCS *.cpp)
FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/gen_tree.h ${CMAKE_SOURCE_DIR}/include/prop_tree.h ${CMAKE_SOURCE_DIR}/include/recon_tree.h)
SET(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
SET(DICTIONARY ${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS} ${HEADERS})
TARGET_LINK_LIBRARIES(${LIBNAME} ${ROOT_LIBRARIES} )

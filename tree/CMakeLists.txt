# Set project name
PROJECT(deuteron_dis_tree)

# find source and header files, set paths, includes...
SET(LIBNAME "deuteron_dis_tree")
FILE(GLOB SRCS *.cpp)
FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/gen_tree.h )
SET(INCLUDES ${ROOT_INCLUDE_DIR} "${CMAKE_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
SET(DICTIONARY ${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS} ${HEADERS})
SET(ROOTSTUFF Gui Gpad Hist Graf Graf3d Tree Rint Postscript Matrix Physics MathCore RIO Net Thread Core Cint Eve Ged RGL TreePlayer EG Geom)
TARGET_LINK_LIBRARIES(${LIBNAME}  "-L${ROOT_LIBRARY_DIR}" ${ROOTSTUFF})

cmake_minimum_required( VERSION 2.6 )
BRAINVISA_PROJECT()

find_package( brainvisa-cmake REQUIRED )

BRAINVISA_DEPENDENCY( DEV DEPENDS "${PROJECT_NAME}" RUN "= ${${PROJECT_NAME}_VERSION}" )

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )
find_package( Doxygen QUIET )

set( AIMS-TIL_INCLUDE_DIRS "${${PROJECT_NAME_UPPER}_BINARY_DIR}/include" ${AIMS-FREE_INCLUDE_DIRS} )
list( REMOVE_DUPLICATES AIMS-TIL_INCLUDE_DIRS )

include( "${AIMS-FREE_USE_FILE}" )

add_subdirectory( src/library )
add_subdirectory( src/commands )

# BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/python"
#                                  ${PROJECT_NAME} )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()


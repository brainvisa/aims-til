cmake_minimum_required( VERSION 2.6 )
BRAINVISA_PROJECT()

find_package( brainvisa-cmake REQUIRED )

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )

include_directories( ${AIMS-FREE_INCLUDE_DIRS} )

add_subdirectory( src/library )
add_subdirectory( src/commands )

# BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/python"
#                                  ${PROJECT_NAME} )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()

BRAINVISA_GET_FILE_LIST_FROM_PRO( cathier_library.pro "HEADERS" _headers "SOURCES" _sources )
set(LIB_HEADERS ${_headers} PARENT_SCOPE)
BRAINVISA_COPY_AND_INSTALL_HEADERS( _headers "" headersTarget )
# BRAINVISA_ADD_MOC_FILES( _sources ${_headers} )
add_library( aimstil SHARED ${_sources} )
add_dependencies( aimstil ${headersTarget} )
set_property( TARGET aimstil PROPERTY VERSION ${${PROJECT_NAME}_VERSION} )
target_link_libraries( aimstil ${_libraries} )

BRAINVISA_INSTALL( TARGETS aimstil
              DESTINATION lib
              COMPONENT ${PROJECT_NAME} )

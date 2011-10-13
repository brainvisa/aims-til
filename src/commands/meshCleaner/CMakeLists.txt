
set( _target "meshCleaner" )
set( _sources "meshCleaner.cc" )
add_executable( ${_target} ${_sources} )
target_link_libraries( ${_target} aimsdata )
BRAINVISA_INSTALL( TARGETS ${_target}
                    DESTINATION bin
                    COMPONENT ${PROJECT_NAME} )

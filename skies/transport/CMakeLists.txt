add_library(skies_transport
	    ${CMAKE_CURRENT_SOURCE_DIR}/elastic.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/iohandler.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/inelastic.cpp    
)
add_dependencies(skies_transport skies_common skies_sampling)

set_target_properties(skies_transport PROPERTIES
                      EXPORT_NAME transport
)

target_link_libraries(skies_transport PUBLIC
	skies_common
	skies_sampling
)

install(TARGETS skies_transport
        DESTINATION lib
	EXPORT skies-targets)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/elastic.h
	${CMAKE_CURRENT_SOURCE_DIR}/iohandler.h
	${CMAKE_CURRENT_SOURCE_DIR}/inelastic.h
	DESTINATION include/skies/transport
)

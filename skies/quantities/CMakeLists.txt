add_library(skies_quantities
	    ${CMAKE_CURRENT_SOURCE_DIR}/dos.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/gmatrix.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/elvelocs.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/eigenvals.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/eigenfreqs.cpp
)
add_dependencies(skies_quantities skies_epw_interface)

set_target_properties(skies_quantities PROPERTIES
		      EXPORT_NAME quantities	
)

target_link_libraries(skies_quantities PUBLIC
	skies_common
	skies_lattices
	skies_sampling
	skies_epw_interface
)

target_include_directories(skies_quantities PUBLIC
	$<INSTALL_INTERFACE:include/skies/quantities>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

install(TARGETS skies_quantities
	DESTINATION lib
	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/dos.h
        ${CMAKE_CURRENT_SOURCE_DIR}/gmatrix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/elvelocs.h
        ${CMAKE_CURRENT_SOURCE_DIR}/eigenvals.h
        ${CMAKE_CURRENT_SOURCE_DIR}/eigenfreqs.h
	DESTINATION include/skies/quantities
)

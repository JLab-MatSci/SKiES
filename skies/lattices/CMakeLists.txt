add_library(skies_lattices
	    ${CMAKE_CURRENT_SOURCE_DIR}/kp_protocol.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/latt_protocol.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/decompose.cpp
)
add_dependencies(skies_lattices skies_common)

set_target_properties(skies_lattices PROPERTIES
                      EXPORT_NAME lattices
)

if(SKIES_ENABLE_MPI)
	target_link_libraries(skies_lattices PUBLIC MPI::MPI_C)
endif()

target_link_libraries(skies_lattices PUBLIC skies_common)

target_include_directories(skies_lattices PUBLIC
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/skies/lattices>
)

install(TARGETS skies_lattices
	DESTINATION lib
	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/kp_protocol.h
        ${CMAKE_CURRENT_SOURCE_DIR}/latt_protocol.h
	DESTINATION include/skies/lattices
)

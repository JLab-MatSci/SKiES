add_library(skies_sampling
			${CMAKE_CURRENT_SOURCE_DIR}/sampling.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/tetrahedra.cpp
)
add_dependencies(skies_sampling skies_common)

set_target_properties(skies_sampling PROPERTIES
                      EXPORT_NAME sampling
)

target_link_libraries(skies_sampling PUBLIC
	skies_common)

target_include_directories(skies_sampling PUBLIC
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/skies/sampling>
)

install(TARGETS skies_sampling
	DESTINATION lib
	EXPORT skies-targets
)

install(FILES
		${CMAKE_CURRENT_SOURCE_DIR}/sampling.h
		${CMAKE_CURRENT_SOURCE_DIR}/tetrahedra.h
	DESTINATION include/skies/sampling
)

add_library(skies_interface INTERFACE)

set_target_properties(skies_interface PROPERTIES
                      EXPORT_NAME interface
)

if (SKIES_ENABLE_MPI)
    add_subdirectory(utils)
    target_link_libraries(skies_interface INTERFACE skies_mpi_wrapper)
endif()

add_subdirectory(common)
target_link_libraries(skies_interface INTERFACE skies_common)

add_subdirectory(lattices)
target_link_libraries(skies_interface INTERFACE skies_lattices)

add_subdirectory(sampling)
target_link_libraries(skies_interface INTERFACE skies_sampling)

add_subdirectory(quantities)
target_link_libraries(skies_interface INTERFACE skies_quantities)

add_subdirectory(interpolation)
target_link_libraries(skies_interface INTERFACE skies_epw_interface)

add_subdirectory(spectral)
target_link_libraries(skies_interface INTERFACE skies_spectral)

add_subdirectory(transport)
target_link_libraries(skies_interface INTERFACE skies_transport)

install(TARGETS skies_interface
        DESTINATION lib
        EXPORT skies-targets
)

install(EXPORT skies-targets
	DESTINATION share/cmake
	NAMESPACE skies::
)

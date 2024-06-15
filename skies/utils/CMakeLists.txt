add_library(skies_mpi_wrapper ${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper.cpp)

set_target_properties(skies_mpi_wrapper PROPERTIES
                      EXPORT_NAME mpi_wrapper
)

target_link_libraries(skies_mpi_wrapper PUBLIC
        MPI::MPI_CXX
)

target_include_directories(skies_mpi_wrapper PUBLIC
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include/skies/utils>
)

install(TARGETS skies_mpi_wrapper
        DESTINATION lib
	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/mpi_wrapper.h
	DESTINATION include/skies/utils
)

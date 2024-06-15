add_library(skies_spectral
	    ${CMAKE_CURRENT_SOURCE_DIR}/lambda.cpp
	    ${CMAKE_CURRENT_SOURCE_DIR}/spec_func.cpp 
)

set_target_properties(skies_spectral PROPERTIES
                      EXPORT_NAME spectral
)

target_link_libraries(skies_spectral PUBLIC skies_quantities)

if (SKIES_ENABLE_MPI)
	target_link_libraries(skies_spectral PUBLIC skies_mpi_wrapper)
endif()

install(TARGETS skies_spectral
        DESTINATION lib
	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/lambda.h
	${CMAKE_CURRENT_SOURCE_DIR}/spec_func.h
	DESTINATION include/skies/spectral
)

add_library(skies_common ${CMAKE_CURRENT_SOURCE_DIR}/alg.cpp)

set_target_properties(skies_common PROPERTIES
		      EXPORT_NAME common
)

target_include_directories(skies_common PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/skies/common>
)

install(TARGETS skies_common
	DESTINATION lib
	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/alg.h
        ${CMAKE_CURRENT_SOURCE_DIR}/units.h
        ${CMAKE_CURRENT_SOURCE_DIR}/ndimarrays.h
	DESTINATION include/skies/common
)

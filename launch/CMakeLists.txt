add_executable(skies
	       ${CMAKE_CURRENT_SOURCE_DIR}/skies.cpp
	       ${CMAKE_CURRENT_SOURCE_DIR}/parser.cpp
	       ${CMAKE_CURRENT_SOURCE_DIR}/commands.cpp
)

if (SKIES_ENABLE_TBB)
	target_link_libraries(skies PRIVATE TBB::tbb)
endif()

target_link_libraries(skies PRIVATE skies_interface)

target_include_directories(skies PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/launch>
)

install(TARGETS skies
    	DESTINATION bin
    	EXPORT skies-targets
)

install(FILES
	${CMAKE_CURRENT_SOURCE_DIR}/parser.h
        ${CMAKE_CURRENT_SOURCE_DIR}/commands.h
		${CMAKE_CURRENT_SOURCE_DIR}/timer.h
        DESTINATION include/launch
)

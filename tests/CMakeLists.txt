set(TEST_SOURCES ndimarrayTest.cpp
                 solveEqBisectTest.cpp
                 customSplitTest.cpp
                 mpiDecompTest.cpp
)

include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG "release-1.11.0"
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/googletest
    UPDATE_COMMAND ""
)
FetchContent_MakeAvailable(googletest)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/googletest/googletest/include)

include(GoogleTest)
add_executable(launch_tests ${TEST_SOURCES})
#set_target_properties(runTests PROPERTIES EXCLUDE_FROM_ALL TRUE)

enable_testing()
target_link_libraries(launch_tests gtest_main skies_interface)

gtest_discover_tests(launch_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

install(TARGETS launch_tests
        DESTINATION bin
)

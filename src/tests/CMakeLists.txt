

# https://rknuus.github.io/blog/2020/04/15/up-and-google-testing-in-minutes
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)


add_executable(tests_lisa 
  ParamBinFileTests.cpp 
  GBSignalTests.cpp)
target_link_libraries(tests_lisa PRIVATE lisa gtest gtest_main)
target_include_directories(tests_lisa PRIVATE include/lisa)
install(TARGETS tests_lisa RUNTIME DESTINATION bin)


target_compile_definitions(tests_lisa 
  PRIVATE 
    $<$<CONFIG:Release>:ENABLE_TIMING_TESTS>
    $<$<CONFIG:RelWithDebInfo>:ENABLE_TIMING_TESTS>)

add_test( tests_lisa tests_lisa )
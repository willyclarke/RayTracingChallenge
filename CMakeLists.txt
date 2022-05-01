project(RayTracingChallenge)

cmake_minimum_required(VERSION 3.0) # setting this is required

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

# if(NOT LINUX) should work, too, if you need that
if(LINUX)
    message(STATUS ">>> Linux")
    # linux stuff here
else()
    message(STATUS ">>> Not Linux")
    # stuff that should happen not on Linux
endif()

##############################################################################
# NOTE: The automatic download of GTest is described on the following web site.
#    -------    https://crascit.com/2015/07/25/cmake-gtest/ ---------
##############################################################################

##############################################################################
# Download and unpack googletest at configure time
##############################################################################
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

##############################################################################
# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
##############################################################################
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

##############################################################################
# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
##############################################################################
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")

##############################################################################
# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
##############################################################################
if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include")
endif()

##############################################################################
# Now simply link your own targets against gtest, gmock,
# etc. as appropriate
##############################################################################

##############################################################################
# The libraries need to be defined first?
##############################################################################
add_subdirectory(common)
add_subdirectory(src/raytrace)

##############################################################################
# Set up the include directories for the targets.
##############################################################################
#target_include_directories(raytrace PUBLIC external)
#target_include_directories(raytrace PUBLIC external/spdlog/include)
target_include_directories(raytrace PUBLIC common/src/main)

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
# The libraries need to be defined first?
add_subdirectory(common)
add_subdirectory(src/raytrace)

# Set up the include directories for the targets.
#target_include_directories(raytrace PUBLIC external)
#target_include_directories(raytrace PUBLIC external/spdlog/include)
target_include_directories(raytrace PUBLIC common/src/main)

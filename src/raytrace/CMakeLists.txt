cmake_minimum_required(VERSION 3.0) # setting this is required
project(raytrace_project)                # this sets the project name
###############################################################################
# Copied from
###############################################################################
# https://github.com/krux02/minimal_cmake_example/blob/master/CMakeLists.txt
###############################################################################
# Also: hints from the following site have been implemented.
# http://voices.canonical.com/jussi.pakkanen/2013/03/26/a-list-of-common-cmake-antipatterns/
###############################################################################

###############################################################################
## file globbing ##############################################################
###############################################################################

# These instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variables
# `sources` and `data`.
file(GLOB_RECURSE sources      src/main/*.cpp src/*.h)
file(GLOB_RECURSE sources_test src/test/*.cpp)
file(GLOB_RECURSE data resources/*)
# You can use set(sources src/main.cpp) etc if you don't want to
# use globbing to find files automatically.

###############################################################################
## Specify debug/release builds ###############################################
###############################################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# NOTE: include the previous setting since it may have been set externally.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -v -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")

###############################################################################
## Specify the c/cpp compiler to use ##########################################
###############################################################################
set(CMAKE_C_COMPILER clang)
#set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER clang++)
#set(CMAKE_CXX_COMPILER g++)

###############################################################################
## definitions ################################################################
###############################################################################
add_definitions(-DHANDMADE_INTERNAL)
###############################################################################
## target definitions #########################################################
###############################################################################

# The data is just added to the executable, because in some IDEs (QtCreator)
# files are invisible when they are not explicitly part of the project.
add_executable(raytrace ${sources} ${data})

# Just for raytrace add some compiler flags.
target_compile_options(raytrace PUBLIC -std=c++1y -Wall -g -Wall -fcolor-diagnostics)
# NOTE: to see more of what is going on you can add "-v"
#target_compile_options(raytrace PUBLIC -v)

# This allows to include files relative to the root of the src directory with a <> pair
target_include_directories(raytrace PUBLIC src/main)

# This copies all resource files in the build directory.
# We need this, because we want to work with paths relative to the executable.
file(COPY ${data} DESTINATION resources)

###############################################################################
## export a compilation database ##############################################
###############################################################################
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

###############################################################################
## dependencies ###############################################################
###############################################################################

# This defines the variables Boost_LIBRARIES that containts all library names
# that we need to link into the program.
#find_package(Boost 1.36.0 COMPONENTS filesystem system REQUIRED)
if (LINUX)
    target_link_libraries(raytrace PUBLIC
        #${Boost_LIBRARIES}
        # here you can add any library dependencies
        pthread
        rt
        raylib
        )
else()
    target_link_libraries(raytrace PUBLIC
        #${Boost_LIBRARIES}
        # here you can add any library dependencies
        pthread
        raylib
        )
endif()

###############################################################################
## testing ####################################################################
###############################################################################

# This is for our testing framework,
# we don't add REQUIRED because it's just for testing.
# People who might want to build the project to use it should not be required
# to install testing dependencies.
find_package(GTest)

if(GTEST_FOUND)
  add_executable(unit_tests ${sources_test} ${sources})

  # This define is added to prevent collision with the main.
  # It might be better solved by not adding the source with the main to the
  # testing target.
  target_compile_definitions(unit_tests PUBLIC UNIT_TESTS)

  # This allows us to use the executable as a link library, and inherit all
  # linker options and library dependencies from it, by simply adding it as dependency.
  set_target_properties(raytrace PROPERTIES ENABLE_EXPORTS on)

  target_link_libraries(unit_tests PUBLIC
    ${GTEST_BOTH_LIBRARIES}
    raytrace
  )

  target_include_directories(unit_tests PUBLIC
    ${GTEST_INCLUDE_DIRS} # doesn't do anything on linux
  )

endif()

###############################################################################
## packaging ##################################################################
###############################################################################

# All install commands get the same destination. this allows us to use paths
# relative to the executable.
install(TARGETS raytrace DESTINATION destination)

# This is basically a repeat of the file copy instruction that copies the
# resources in the build directory, but here we tell cmake that we want it
# in the package.
#install(DIRECTORY resources DESTINATION raytrace_destination)

# Now comes everything we need, to create a package
# there are a lot more variables you can set, and some
# you need to set for some package types, but we want to
# be minimal here.
set(CPACK_PACKAGE_NAME "Myraytrace")
set(CPACK_PACKAGE_VERSION "1.0.0")

# We don't want to split our program up into several incomplete pieces.
set(CPACK_MONOLITHIC_INSTALL 1)

# This must be last
include(CPack)


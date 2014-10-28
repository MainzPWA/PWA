cmake_minimum_required (VERSION 2.8)
project(PWA)
set(PROJ "PWA")

# check for in-source build, forbid it!
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR AND NOT MSVC_IDE)
  message(FATAL_ERROR "\nIn-source build attempt detected!\n"
    "Please create a new directory (e.g. build) and run `cmake ..`.\n"
    "Also don't forget to delete the created CMakeCache.txt and CMakeFiles dir"
    )
endif()

message(STATUS "*** Build Type:      " ${CMAKE_BUILD_TYPE})
message(STATUS "*** Compiler Flags:  " ${DEFAULT_COMPILE_FLAGS})
message(STATUS "*** Install libs to: " ${LIBRARY_OUTPUT_PATH})
message(STATUS "*** Install bin to:  " ${EXECUTABLE_OUTPUT_PATH})

# use some sane default directories
if(NOT DEFINED EXECUTABLE_OUTPUT_PATH)
       set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
endif()

# we check for empty string here, since the variable is indeed defined to an empty string
if(NOT CMAKE_BUILD_TYPE)
  # this also reflects this default in the GUI
        SET(CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# really no optimization in debug mode
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O0")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -s")
endif()

# require a fairly recent ROOT version
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(ROOT "5.28" REQUIRED)
link_directories(${ROOT_LIBRARY_DIR})

# define command for the IncreaseBuild executable (increases build number in Build.h)
add_executable(IncreaseBuild ${CMAKE_SOURCE_DIR}/cmake/IncreaseBuild.c)
# define  target for header file carrying automatic counted build number
add_custom_target(Build COMMAND IncreaseBuild ${CMAKE_SOURCE_DIR}/src/Build.h)

# construct the executable from source files
include_directories("src" ${ROOT_INCLUDE_DIR})
aux_source_directory("src" SRC_LIST)
add_executable(${PROJ} ${SRC_LIST})
# make project's main executable dependend on build number header file
add_dependencies(${PROJ} Build)
target_link_libraries(${PROJ} ${ROOT_LIBRARIES} "-lMinuit")

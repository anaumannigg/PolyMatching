cmake_minimum_required(VERSION 3.3)

project(PolygonMatching)

set(CMAKE_CXX_STANDARD 17)

option(CXX "enable C++ compilation" ON)

if(CXX)
  enable_language(CXX)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI REQUIRED)

find_package(CGAL REQUIRED)

find_package(Threads REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# list source files here
set(sources
        src/main.cpp
  src/shapefile_io_operations.cpp
  src/localization.cpp
  src/graph_computations.cpp
  src/solution.cpp
  src/linear_program.cpp
        include/cgal_includes.h
  include/shapefile_io_operations.h
  include/localization.h
  include/graph_computations.h
  include/solution.h
  include/linear_program.h
        include/command_line_parser.h
        src/command_line_parser.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})
target_link_libraries(${CMAKE_PROJECT_NAME} CGAL::CGAL)
add_subdirectory(lib/shapelib)
target_link_libraries(${CMAKE_PROJECT_NAME} shp)
target_link_libraries(${CMAKE_PROJECT_NAME} Threads::Threads)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()

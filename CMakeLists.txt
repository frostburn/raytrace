CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
PROJECT(raytrace CXX)

SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
SET(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/include")
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})

ADD_SUBDIRECTORY(src)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)

ADD_EXECUTABLE(render main.cpp)
TARGET_LINK_LIBRARIES(render raytrace)

ADD_EXECUTABLE(render_cartesian main_cartesian.cpp)
TARGET_LINK_LIBRARIES(render_cartesian raytrace)

ADD_EXECUTABLE(debug_sphere debug_sphere.cpp)
TARGET_LINK_LIBRARIES(debug_sphere raytrace)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra --std c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)
ADD_LIBRARY(
  raytrace
  quaternion.cpp
  ray_traceable.cpp
  gradient_traceable.cpp
  approach_traceable.cpp
  sphere.cpp
  clifford_torus.cpp
  trace.cpp
  input_parser.cpp
  hyper_torus.cpp
)
INSTALL(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include/raytrace
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
)
INSTALL(
  TARGETS raytrace
  ARCHIVE
  DESTINATION lib
)

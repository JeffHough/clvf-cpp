project(clvf)
cmake_minimum_required(VERSION 3.16)

# require Eigen3 to be installed:
find_package(Eigen3 REQUIRED)

# create the spacecraft/clvf library:
add_library(
  clvf
  STATIC
  src/clvf.cc
  src/spacecraft.cc
)

# link to eigen:
target_link_libraries(clvf Eigen3::Eigen)

# the include directories:
target_include_directories(clvf PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
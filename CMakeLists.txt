cmake_minimum_required(VERSION 3.14)
project(cosmology_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W4 /utf-8)  # sources en UTF-8 (accents)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

add_library(cosmology STATIC Cosmology.cpp SvgPlot.cpp Report.cpp)
target_include_directories(cosmology PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(cosmology_simulation main.cpp)
target_link_libraries(cosmology_simulation PRIVATE cosmology)

enable_testing()
add_executable(test_cosmology tests/test_cosmology.cpp)
target_link_libraries(test_cosmology PRIVATE cosmology)
add_test(NAME physique COMMAND test_cosmology)

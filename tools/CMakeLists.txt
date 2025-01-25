CMAKE_MINIMUM_REQUIRED ( VERSION 3.10.0 )
include_directories(.)
include_directories(../src/core)

add_subdirectory(ext/glog)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(brdf simbrdf.cpp hairsimbrdf.cpp spectrum.h spectrum.cpp)

target_link_libraries(brdf glog::glog)
target_link_libraries(brdf ${OpenCV_LIBS})
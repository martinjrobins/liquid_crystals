cmake_minimum_required(VERSION 2.8)
find_package(Boost 1.50.0 COMPONENTS filesystem system date_time python REQUIRED)

find_package(PythonLibs REQUIRED)

find_package(VTK REQUIRED)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x -fopenmp -D_GLIBCXX_PARALLEL")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

include_directories(${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR}  ${PYTHON_INCLUDE_DIRS} src Aboria/src Aboria/eigen3)

add_library (ParticleSimulation src/Python.cpp)
TARGET_LINK_LIBRARIES(ParticleSimulation ${PYTHON_LIBRARIES} ${VTK_LIBRARIES})
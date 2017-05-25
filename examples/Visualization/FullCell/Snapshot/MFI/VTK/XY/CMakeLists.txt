cmake_minimum_required(VERSION 2.6)
project(Test)
set(VTK_DIR "/usr/local/share/cmake-3.0/Modules")
find_package(VTK REQUIRED)
find_package(Qt5Widgets)
include(${VTK_USE_FILE})
add_executable(simulate main.cxx)
target_link_libraries(simulate ${VTK_LIBRARIES})


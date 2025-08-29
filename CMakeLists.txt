cmake_minimum_required(VERSION 3.16)
project(DysonSwarm CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)

add_executable(DysonSwarm main.cpp)

target_link_libraries(DysonSwarm
  PRIVATE
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
    imgui::imgui
)

# Activa docking/viewports en el código
target_compile_definitions(DysonSwarm PRIVATE IMGUI_ENABLE_DOCKING IMGUI_ENABLE_VIEWPORTS)

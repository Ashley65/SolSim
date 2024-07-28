cmake_minimum_required(VERSION 3.28)
project(SolSim CXX)

set(CMAKE_CXX_STANDARD 26)


# Adding python3 to the project
find_package(Python3 COMPONENTS Interpreter Development)

set(PHYSICS_ENGINE_INCLUDE_DIR
        include/physicsEngine/simulator.h
        include/physicsEngine/CelestialBody.h
        include/physicsEngine/RKF45.h)

set(PHYSICS_ENGINE_SOURCE_DIR

        src/physicsEngine/CelestialBody.cpp
        src/physicsEngine/RKF45.cpp
        )


add_executable(SolSim main.cpp
        src/app.cpp
        src/app.h
        ${PHYSICS_ENGINE_INCLUDE_DIR}
        ${PHYSICS_ENGINE_SOURCE_DIR})





find_package(OpenGL REQUIRED)
target_link_libraries(SolSim PRIVATE OpenGL::GL)

find_package(glm CONFIG REQUIRED)
target_link_libraries(SolSim PRIVATE glm::glm)

find_package(GLEW REQUIRED)
target_link_libraries(SolSim PRIVATE GLEW::GLEW)


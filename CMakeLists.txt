cmake_minimum_required(VERSION 3.12)
project(untitled)
# set(CMAKE_BUILD_TYPE debug)
set(CMAKE_BUILD_TYPE release)
set(CPP_LIBRARIES_PATH "C:/Users/xd15zhn/Documents/cpplibraries")
list(APPEND CMAKE_PREFIX_PATH ${CPP_LIBRARIES_PATH})

set(PROJECT_SOURCES
    ${PROJECT_SOURCE_DIR}/main.cpp
    ${PROJECT_SOURCE_DIR}/src/detector.cpp
    ${PROJECT_SOURCE_DIR}/src/zhnoptim.cpp
    ${PROJECT_SOURCE_DIR}/src/opt_de.cpp
    ${PROJECT_SOURCE_DIR}/src/opt_ps.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})

find_package(simucpp REQUIRED)
message(STATUS "simucpp_LIBS: ${simucpp_LIBS}")
message(STATUS "simucpp_VERSION: ${simucpp_VERSION}")
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CPP_LIBRARIES_PATH}/include
    ${PROJECT_SOURCE_DIR}/src
    ${simucpp_INCLUDE_DIRS}
)
target_link_libraries(${PROJECT_NAME} PUBLIC
    ${simucpp_LIBS}
)

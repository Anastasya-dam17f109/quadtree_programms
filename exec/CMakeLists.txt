cmake_minimum_required(VERSION 3.9)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
project(exec)

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

find_package(OpenMP  REQUIRED)
find_package(jsoncpp CONFIG REQUIRED)
# get_target_property(JSON_INC_PATH jsoncpp_lib  INTERFACE_INCLUDE_DIRECTORIES)
# include_directories({JSON_INC_PATH})


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME}
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(${PROJECT_NAME} PUBLIC 
Handlers 
OpenMP::OpenMP_CXX 
jsoncpp_lib 
jsoncpp_object 
JsonCpp::JsonCpp
)

cmake_minimum_required(VERSION 3.23)
project(Pathtracing)

include(FetchContent)

file(GLOB sources src/*.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

add_executable(main ${sources})
target_compile_options(main PUBLIC -std=c++17 -Wall -Wextra -Wconversion -O3)
include_directories(include)

FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(json)

target_link_libraries(main PRIVATE nlohmann_json::nlohmann_json)

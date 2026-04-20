cmake_minimum_required(VERSION 3.23)
project(Pathtracing)

include(FetchContent)

file(GLOB sources src/*.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

add_executable(pathtracer ${sources} app/pathtracer.cpp)
target_compile_options(pathtracer PUBLIC -std=c++17 -Wall -Wextra -Wconversion -Werror -O3)

add_executable(denoiser ${sources} app/denoiser.cpp)
target_compile_options(denoiser PUBLIC -std=c++17 -Wall -Wextra -Wconversion -Werror -O3)


include_directories(include)

FetchContent_Declare(
  json
  URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_MakeAvailable(json)

target_link_libraries(pathtracer PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(denoiser PRIVATE nlohmann_json::nlohmann_json)

cmake_minimum_required(VERSION 3.14)
project(numeric
    VERSION 1.0.0
    DESCRIPTION "Arbitrary-precision integer and floating-point arithmetic in C++"
    LANGUAGES CXX
)

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Options ───────────────────────────────────────────────────────────────────
option(NUMERIC_BUILD_TESTS "Build unit tests" ON)
option(NUMERIC_BUILD_MAIN "Build the demo/benchmark main" ON)

# ── Testing (must be in root for add_test in subdirs) ──────────────────────
if(NUMERIC_BUILD_TESTS)
    enable_testing()
endif()

# ── Threading ─────────────────────────────────────────────────────────────────
find_package(Threads REQUIRED)

# ── Interface library (header-only) ───────────────────────────────────────────
add_library(numeric-lib INTERFACE)
target_include_directories(numeric-lib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(numeric-lib INTERFACE Threads::Threads)

# ── Demo / Benchmark executable ───────────────────────────────────────────────
if(NUMERIC_BUILD_MAIN)
    add_executable(numeric_demo src/main.cpp)
    target_link_libraries(numeric_demo PRIVATE numeric-lib)
    target_compile_features(numeric_demo PRIVATE cxx_std_23)
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────
if(NUMERIC_BUILD_TESTS)
    add_subdirectory(tests)
endif()

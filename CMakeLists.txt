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

include(CheckCXXSourceCompiles)

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

# ── AVX/AVX2 compiler support detection ──────────────────────────────────────
set(CMAKE_REQUIRED_FLAGS "-mavx")
check_cxx_source_compiles("
    #include <immintrin.h>
    int main() {
        __m128i a = _mm_setzero_si128();
        __m128i b = _mm_add_epi64(a, a);
        (void)b;
        return 0;
    }
" NN_HAVE_AVX_SUPPORT)
set(CMAKE_REQUIRED_FLAGS "")

set(CMAKE_REQUIRED_FLAGS "-mavx2")
check_cxx_source_compiles("
    #include <immintrin.h>
    int main() {
        __m256i a = _mm256_setzero_si256();
        __m256i b = _mm256_add_epi64(a, a);
        (void)b;
        return 0;
    }
" NN_HAVE_AVX2_SUPPORT)
set(CMAKE_REQUIRED_FLAGS "")

# ── Demo / Benchmark executable ───────────────────────────────────────────────
if(NUMERIC_BUILD_MAIN)
    add_executable(numeric_demo src/main.cpp)
    target_link_libraries(numeric_demo PRIVATE numeric-lib)
    target_compile_features(numeric_demo PRIVATE cxx_std_23)

    # ── AVX variant ────────────────────────────────────────────────────────
    if(NN_HAVE_AVX_SUPPORT)
        add_executable(numeric_demo_avx src/main.cpp)
        target_compile_options(numeric_demo_avx PRIVATE -mavx -DNN_HAVE_AVX)
        target_link_libraries(numeric_demo_avx PRIVATE numeric-lib)
        target_compile_features(numeric_demo_avx PRIVATE cxx_std_23)
        message(STATUS "AVX support detected: building numeric_demo_avx")
    else()
        message(STATUS "AVX not supported by compiler: skipping numeric_demo_avx")
    endif()

    # ── AVX2 variant ───────────────────────────────────────────────────────
    if(NN_HAVE_AVX2_SUPPORT)
        add_executable(numeric_demo_avx2 src/main.cpp)
        target_compile_options(numeric_demo_avx2 PRIVATE -mavx2 -DNN_HAVE_AVX2)
        target_link_libraries(numeric_demo_avx2 PRIVATE numeric-lib)
        target_compile_features(numeric_demo_avx2 PRIVATE cxx_std_23)
        message(STATUS "AVX2 support detected: building numeric_demo_avx2")
    else()
        message(STATUS "AVX2 not supported by compiler: skipping numeric_demo_avx2")
    endif()
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────
if(NUMERIC_BUILD_TESTS)
    add_subdirectory(tests)
endif()

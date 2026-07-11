# ── Integer test ──────────────────────────────────────────────────────────────
add_executable(test_integer test_integer.cpp)
target_link_libraries(test_integer PRIVATE numeric-lib)
target_compile_features(test_integer PRIVATE cxx_std_17)
add_test(NAME test_integer COMMAND test_integer)

# ── Numeric test ──────────────────────────────────────────────────────────────
add_executable(test_numeric test_numeric.cpp)
target_link_libraries(test_numeric PRIVATE numeric-lib)
target_compile_features(test_numeric PRIVATE cxx_std_17)
add_test(NAME test_numeric COMMAND test_numeric)

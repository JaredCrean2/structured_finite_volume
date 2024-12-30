enable_testing()
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})

# Test utils
file(GLOB utils_src_files utils/*.cc)
add_library(test_utils "${utils_src_files}")
target_include_directories(test_utils PUBLIC
                           "${CMAKE_CURRENT_SOURCE_DIR}/utils"
)

target_link_libraries(test_utils
                     mesh
                     ${GTEST_BOTH_LIBRARIES}
                   )
 
# Unit tests
file(GLOB unit_src_files unit/*.cc)

add_executable(unit_tests ${unit_src_files})
target_include_directories(unit_tests PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/src"
)


target_link_libraries(unit_tests
                     ${GTEST_BOTH_LIBRARIES}
                     pthread 
                     mesh
                     test_utils)


# integration tests
file(GLOB integration_src_files integration/*.cc)

add_executable(integration_tests ${integration_src_files})
target_include_directories(unit_tests PUBLIC
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/src"
)


target_link_libraries(integration_tests
                     ${GTEST_BOTH_LIBRARIES}
                     pthread 
                     mesh
                     test_utils)                  

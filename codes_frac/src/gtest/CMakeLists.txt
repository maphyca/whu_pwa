set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/gtest/)
add_executable(test_dataobj data_obj_unittest.cc ${PROJECT_SOURCE_DIR}/src/cpp/data_obj_new.cc)
target_link_libraries(test_dataobj gtest gtest_main)

aux_source_directory(. PWA_SOURCES)
add_executable(test_m ${PWA_SOURCES})
target_link_libraries(test_m ROOT::Minuit2 test_cuda)

# define the test executable:
add_executable(test test.cc)

# link it to the clvf library:
target_link_libraries(test PRIVATE clvf)
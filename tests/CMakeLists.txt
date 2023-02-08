# define the test executable:
add_executable(test test.cc)
add_executable(consistent_transformation_test consistent_transformations.cc)

# link it to the clvf library:
target_link_libraries(test PRIVATE clvf)
target_link_libraries(consistent_transformation_test PRIVATE clvf)

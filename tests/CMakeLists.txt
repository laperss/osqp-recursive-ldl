# Add subdirectories of all the tests
# add_subdirectory(basic_qp)
# add_subdirectory(basic_qp2)
# add_subdirectory(lin_alg)
# add_subdirectory(non_cvx)
# add_subdirectory(primal_dual_infeasibility)
# add_subdirectory(primal_infeasibility)
# add_subdirectory(solve_linsys)
# add_subdirectory(unconstrained)
# add_subdirectory(update_matrices)

# set(test_headers ${headers} PARENT_SCOPE)
# set(codegen_test_headers ${codegen_headers} PARENT_SCOPE)





# Google Test
add_executable(
    unit_tests
    example_add.cpp
    example_subtract.cpp
    )

target_link_libraries(
    unit_tests
    gtest_main
    example
    )

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/unit_tests
  )

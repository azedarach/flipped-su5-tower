set(cSMHdCKM_TEST_TARGETS
  test_cSMHdCKM_Kappa_running
)

set(cSMHdCKMRHNEFT_TEST_TARGETS
  test_cSMHdCKMRHNEFT_boundary_conditions
  test_cSMHdCKMRHNEFT_low_scale_constraint
)

foreach(test_target ${cSMHdCKM_TEST_TARGETS})
  add_executable(
    ${test_target}
    ${test_target}.cpp
    )

  target_include_directories(
    ${test_target}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

  target_link_libraries(
    ${test_target}
    FlexibleSUSY::cSMHdCKM
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )

  if(ENABLE_COMPILER_WARNINGS)
    target_compile_options(${test_target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
  endif(ENABLE_COMPILER_WARNINGS)

  add_test(NAME ${test_target} COMMAND ${test_target})
endforeach(test_target)

foreach(test_target ${cSMHdCKMRHNEFT_TEST_TARGETS})
  add_executable(
    ${test_target}
    ${test_target}.cpp
    )

  target_include_directories(
    ${test_target}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

  target_link_libraries(
    ${test_target}
    FlexibleSUSY::cSMHdCKM
    FlexibleSUSY::cSMHdCKMRHN
    FlexibleSUSY::cSMHdCKMRHNEFT
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )

  if(ENABLE_COMPILER_WARNINGS)
    target_compile_options(${test_target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
  endif(ENABLE_COMPILER_WARNINGS)

  add_test(NAME ${test_target} COMMAND ${test_target})
endforeach(test_target)

set(ALL_TEST_TARGETS
  ${cSMHdCKM_TEST_TARGETS}
  ${cSMHdCKMRHNEFT_TEST_TARGETS}
)


add_custom_target(check COMMAND "ctest")
add_dependencies(check ${ALL_TEST_TARGETS})

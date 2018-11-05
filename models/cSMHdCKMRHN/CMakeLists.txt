set(cSMHdCKMRHN_SRC
  cSMHdCKMRHN_a_muon.cpp
  cSMHdCKMRHN_edm.cpp
  cSMHdCKMRHN_effective_couplings.cpp
  cSMHdCKMRHN_info.cpp
  cSMHdCKMRHN_input_parameters.cpp
  cSMHdCKMRHN_mass_eigenstates.cpp
  cSMHdCKMRHN_observables.cpp
  cSMHdCKMRHN_physical.cpp
  cSMHdCKMRHN_slha_io.cpp
  cSMHdCKMRHN_soft_beta_mu2.cpp
  cSMHdCKMRHN_soft_beta_Mv.cpp
  cSMHdCKMRHN_soft_beta_v.cpp
  cSMHdCKMRHN_soft_parameters.cpp
  cSMHdCKMRHN_susy_beta_g1.cpp
  cSMHdCKMRHN_susy_beta_g2.cpp
  cSMHdCKMRHN_susy_beta_g3.cpp
  cSMHdCKMRHN_susy_beta_Lambdax.cpp
  cSMHdCKMRHN_susy_beta_Yd.cpp
  cSMHdCKMRHN_susy_beta_Ye.cpp
  cSMHdCKMRHN_susy_beta_Yu.cpp
  cSMHdCKMRHN_susy_beta_Yv.cpp
  cSMHdCKMRHN_susy_parameters.cpp
  cSMHdCKMRHN_utilities.cpp
  cSMHdCKMRHN_weinberg_angle.cpp
)

if(ENABLE_TWO_SCALE_SOLVER)
  set(cSMHdCKMRHN_SRC
    ${cSMHdCKMRHN_SRC}
    cSMHdCKMRHN_two_scale_convergence_tester.cpp
    cSMHdCKMRHN_two_scale_ewsb_solver.cpp
    cSMHdCKMRHN_two_scale_high_scale_constraint.cpp
    cSMHdCKMRHN_two_scale_initial_guesser.cpp
    cSMHdCKMRHN_two_scale_low_scale_constraint.cpp
    cSMHdCKMRHN_two_scale_model.cpp
    cSMHdCKMRHN_two_scale_spectrum_generator.cpp
    cSMHdCKMRHN_two_scale_susy_scale_constraint.cpp
  )
endif(ENABLE_TWO_SCALE_SOLVER)

if(ENABLE_SEMI_ANALYTIC_SOLVER)

endif(ENABLE_SEMI_ANALYTIC_SOLVER)

add_library(cSMHdCKMRHN ${cSMHdCKMRHN_SRC})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(cSMHdCKMRHN PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(cSMHdCKMRHN
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cSMHdCKMRHN
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
)

if(ENABLE_LIBRARYLINK)
  set(LLcSMHdCKMRHN_SRC
    cSMHdCKMRHN_librarylink.cpp
    )

  add_library(LLcSMHdCKMRHN ${LLcSMHdCKMRHN_SRC})

  if(ENABLE_COMPILER_WARNINGS)
    target_compile_options(LLcSMHdCKMRHN PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
  endif(ENABLE_COMPILER_WARNINGS)

  target_include_directories(LLcSMHdCKMRHN
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Mathematica_MathLink_INCLUDE_DIR}
    ${Mathematica_WolframLibrary_INCLUDE_DIR}
    )

  target_link_libraries(LLcSMHdCKMRHN
    PUBLIC
    Boost::boost
    Eigen3::Eigen
    FlexibleSUSY::config
    FlexibleSUSY::fflite
    FlexibleSUSY::slhaea
    FlexibleSUSY::flexisusy
    FlexibleSUSY::model_specific::SM
    GSL::gsl
    LAPACK::LAPACK
    LoopTools::LoopTools
    ${Mathematica_MathLink_LIBRARY}
    SQLite3::SQLite3
    Threads::Threads
    )

endif(ENABLE_LIBRARYLINK)

add_executable(run_cSMHdCKMRHN
  run_cSMHdCKMRHN.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_cSMHdCKMRHN PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_cSMHdCKMRHN
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKMRHN
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
)

add_executable(run_cmd_line_cSMHdCKMRHN
  run_cmd_line_cSMHdCKMRHN.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_cmd_line_cSMHdCKMRHN PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_cmd_line_cSMHdCKMRHN
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKMRHN
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
)

add_executable(scan_cSMHdCKMRHN
  scan_cSMHdCKMRHN.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(scan_cSMHdCKMRHN PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(scan_cSMHdCKMRHN
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKMRHN
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
)

add_library(FlexibleSUSY::cSMHdCKMRHN ALIAS cSMHdCKMRHN)
if(ENABLE_LIBRARYLINK)
  add_library(FlexibleSUSY::LLcSMHdCKMRHN ALIAS LLcSMHdCKMRHN)
endif(ENABLE_LIBRARYLINK)

add_executable(FlexibleSUSY::run_cSMHdCKMRHN ALIAS run_cSMHdCKMRHN)
add_executable(FlexibleSUSY::run_cmd_line_cSMHdCKMRHN ALIAS run_cmd_line_cSMHdCKMRHN)
add_executable(FlexibleSUSY::scan_cSMHdCKMRHN ALIAS scan_cSMHdCKMRHN)

set(cSMHdCKM_SRC
  cSMHdCKM_a_muon.cpp
  cSMHdCKM_edm.cpp
  cSMHdCKM_effective_couplings.cpp
  cSMHdCKM_info.cpp
  cSMHdCKM_input_parameters.cpp
  cSMHdCKM_mass_eigenstates.cpp
  cSMHdCKM_observables.cpp
  cSMHdCKM_physical.cpp
  cSMHdCKM_slha_io.cpp
  cSMHdCKM_soft_beta_mu2.cpp
  cSMHdCKM_soft_beta_v.cpp
  cSMHdCKM_soft_parameters.cpp
  cSMHdCKM_susy_beta_g1.cpp
  cSMHdCKM_susy_beta_g2.cpp
  cSMHdCKM_susy_beta_g3.cpp
  cSMHdCKM_susy_beta_Lambdax.cpp
  cSMHdCKM_susy_beta_Yd.cpp
  cSMHdCKM_susy_beta_Ye.cpp
  cSMHdCKM_susy_beta_Yu.cpp
  cSMHdCKM_susy_parameters.cpp
  cSMHdCKM_utilities.cpp
  cSMHdCKM_weinberg_angle.cpp
)

if(ENABLE_TWO_SCALE_SOLVER)
  set(cSMHdCKM_SRC
    ${cSMHdCKM_SRC}
    cSMHdCKM_two_scale_convergence_tester.cpp
    cSMHdCKM_two_scale_ewsb_solver.cpp
    cSMHdCKM_two_scale_high_scale_constraint.cpp
    cSMHdCKM_two_scale_initial_guesser.cpp
    cSMHdCKM_two_scale_low_scale_constraint.cpp
    cSMHdCKM_two_scale_model.cpp
    cSMHdCKM_two_scale_spectrum_generator.cpp
    cSMHdCKM_two_scale_susy_scale_constraint.cpp
  )
endif(ENABLE_TWO_SCALE_SOLVER)

if(ENABLE_SEMI_ANALYTIC_SOLVER)

endif(ENABLE_SEMI_ANALYTIC_SOLVER)

add_library(cSMHdCKM ${cSMHdCKM_SRC})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(cSMHdCKM PRIVATE ${CXX_WARNING_FLAGS})
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(cSMHdCKM
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cSMHdCKM
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
  ${SQLITE3_LIBRARY}
  Threads::Threads
)

if(ENABLE_LIBRARYLINK)
  set(LLcSMHdCKM_SRC
    cSMHdCKM_librarylink.cpp
    )

  add_library(LLcSMHdCKM ${LLcSMHdCKM_SRC})

  if(ENABLE_COMPILER_WARNINGS)
    target_compile_options(LLcSMHdCKM PRIVATE ${CXX_WARNING_FLAGS})
  endif(ENABLE_COMPILER_WARNINGS)

  target_include_directories(LLcSMHdCKM
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Mathematica_MathLink_INCLUDE_DIR}
    ${Mathematica_WolframLibrary_INCLUDE_DIR}
    )

  target_link_libraries(LLcSMHdCKM
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
    ${SQLITE3_LIBRARY}
    Threads::Threads
    )

endif(ENABLE_LIBRARYLINK)

add_executable(run_cSMHdCKM
  run_cSMHdCKM.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_cSMHdCKM PRIVATE ${CXX_WARNING_FLAGS})
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_cSMHdCKM
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKM
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  ${SQLITE3_LIBRARY}
  Threads::Threads
)

add_executable(run_cmd_line_cSMHdCKM
  run_cmd_line_cSMHdCKM.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_cmd_line_cSMHdCKM PRIVATE ${CXX_WARNING_FLAGS})
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_cmd_line_cSMHdCKM
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKM
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  ${SQLITE3_LIBRARY}
  Threads::Threads
)

add_executable(scan_cSMHdCKM
  scan_cSMHdCKM.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(scan_cSMHdCKM PRIVATE ${CXX_WARNING_FLAGS})
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(scan_cSMHdCKM
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKM
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  ${SQLITE3_LIBRARY}
  Threads::Threads
)

add_library(FlexibleSUSY::cSMHdCKM ALIAS cSMHdCKM)
if(ENABLE_LIBRARYLINK)
  add_library(FlexibleSUSY::LLcSMHdCKM ALIAS LLcSMHdCKM)
endif(ENABLE_LIBRARYLINK)

add_executable(FlexibleSUSY::run_cSMHdCKM ALIAS run_cSMHdCKM)
add_executable(FlexibleSUSY::run_cmd_line_cSMHdCKM ALIAS run_cmd_line_cSMHdCKM)
add_executable(FlexibleSUSY::scan_cSMHdCKM ALIAS scan_cSMHdCKM)

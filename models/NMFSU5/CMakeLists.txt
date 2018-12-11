set(NMFSU5_src
  NMFSU5_info.cpp
  NMFSU5_input_parameters.cpp
  NMFSU5_mass_eigenstates.cpp
  NMFSU5_physical.cpp
  NMFSU5_soft_beta_m10sq.cpp
  NMFSU5_soft_beta_m12sq.cpp
  NMFSU5_soft_beta_m5Prsq.cpp
  NMFSU5_soft_beta_m5sq.cpp
  NMFSU5_soft_beta_mu.cpp
  NMFSU5_soft_beta_muPr.cpp
  NMFSU5_soft_beta_v.cpp
  NMFSU5_soft_beta_VG.cpp
  NMFSU5_soft_beta_vPr.cpp
  NMFSU5_soft_parameters.cpp
  NMFSU5_susy_beta_eta1.cpp
  NMFSU5_susy_beta_eta2.cpp
  NMFSU5_susy_beta_eta3.cpp
  NMFSU5_susy_beta_g5.cpp
  NMFSU5_susy_beta_gX.cpp
  NMFSU5_susy_beta_lam1.cpp
  NMFSU5_susy_beta_lam2.cpp
  NMFSU5_susy_beta_lam3.cpp
  NMFSU5_susy_beta_lam3t.cpp
  NMFSU5_susy_beta_lam4.cpp
  NMFSU5_susy_beta_lam4t.cpp
  NMFSU5_susy_beta_lam5.cpp
  NMFSU5_susy_beta_lam5t.cpp
  NMFSU5_susy_beta_lam6.cpp
  NMFSU5_susy_beta_lam6t.cpp
  NMFSU5_susy_beta_lam7.cpp
  NMFSU5_susy_beta_lam8.cpp
  NMFSU5_susy_beta_Y10.cpp
  NMFSU5_susy_beta_Y10Pr.cpp
  NMFSU5_susy_beta_Y1.cpp
  NMFSU5_susy_beta_Y1Pr.cpp
  NMFSU5_susy_beta_Y5b.cpp
  NMFSU5_susy_beta_Y5bPr.cpp
  NMFSU5_susy_parameters.cpp
)

if(ENABLE_TWO_SCALE_SOLVER)
  set(cSMHdCKMRHN_SRC
    ${cSMHdCKMRHN_SRC}
    cSMHdCKMRHN_two_scale_ewsb_solver.cpp
  )
endif(ENABLE_TWO_SCALE_SOLVER)

add_library(NMFSU5 ${NMFSU5_src})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(NMFSU5 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(NMFSU5
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(NMFSU5
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

add_library(FlexibleSUSY::NMFSU5 ALIAS NMFSU5)

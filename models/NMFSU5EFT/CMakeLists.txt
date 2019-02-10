set(NMFSU5EFT_SRC
  NMFSU5EFT_input_parameters.cpp
  NMFSU5EFT_matching.cpp
  NMFSU5EFT_slha_io.cpp
  )

if(ENABLE_TWO_SCALE_SOLVER)
  set(NMFSU5EFT_SRC
    ${NMFSU5EFT_SRC}
    NMFSU5EFT_two_scale_high_scale_constraint.cpp
    NMFSU5EFT_two_scale_initial_guesser.cpp
    NMFSU5EFT_two_scale_low_scale_constraint.cpp
    NMFSU5EFT_two_scale_matching.cpp
    NMFSU5EFT_two_scale_spectrum_generator.cpp
    NMFSU5EFT_two_scale_susy_scale_constraint.cpp
    )
endif(ENABLE_TWO_SCALE_SOLVER)

if(ENABLE_SEMI_ANALYTIC_SOLVER)

endif(ENABLE_SEMI_ANALYTIC_SOLVER)

add_library(cSMHdCKMRHNEFT ${cSMHdCKMRHNEFT_SRC})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(NMFSU5EFT PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(NMFSU5EFT
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(NMFSU5EFT
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  FlexibleSUSY::cSMHdCKM
  FlexibleSUSY::cSMHdCKMRHN
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
  )

add_executable(run_NMFSU5EFT
  run_NMFSU5EFT.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_NMFSU5EFT PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_NMFSU5EFT
  PUBLIC
  Boost::boost
  Eigen3::Eigen
  FlexibleSUSY::config
  FlexibleSUSY::fflite
  FlexibleSUSY::slhaea
  FlexibleSUSY::flexisusy
  FlexibleSUSY::model_specific::SM
  cSMHdCKMRHNEFT
  GSL::gsl
  LAPACK::LAPACK
  LoopTools::LoopTools
  SQLite3::SQLite3
  Threads::Threads
)

add_library(FlexibleSUSY::NMFSU5EFT ALIAS NMFSU5EFT)

add_executable(FlexibleSUSY::run_NMFSU5EFT ALIAS run_NMFSU5EFT)

set(cSMHdCKMRHNEFT_SRC
  cSMHdCKMRHNEFT_matching.cpp
  cSMHdCKMRHNEFT_slha_io.cpp
  )

if(ENABLE_TWO_SCALE_SOLVER)
  set(cSMHdCKMRHNEFT_SRC
    ${cSMHdCKMRHNEFT_SRC}
    cSMHdCKMRHNEFT_two_scale_initial_guesser.cpp
    cSMHdCKMRHNEFT_two_scale_matching.cpp
    cSMHdCKMRHNEFT_two_scale_spectrum_generator.cpp
    )
endif(ENABLE_TWO_SCALE_SOLVER)

if(ENABLE_SEMI_ANALYTIC_SOLVER)

endif(ENABLE_SEMI_ANALYTIC_SOLVER)

add_library(cSMHdCKMRHNEFT ${cSMHdCKMRHNEFT_SRC})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(cSMHdCKMRHNEFT PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(cSMHdCKMRHNEFT
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cSMHdCKMRHNEFT
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

add_executable(run_cSMHdCKMRHNEFT
  run_cSMHdCKMRHNEFT.cpp
)

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(run_cSMHdCKMRHNEFT PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_link_libraries(run_cSMHdCKMRHNEFT
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

add_library(FlexibleSUSY::cSMHdCKMRHNEFT ALIAS cSMHdCKMRHNEFT)

add_executable(FlexibleSUSY::run_cSMHdCKMRHNEFT ALIAS run_cSMHdCKMRHNEFT)

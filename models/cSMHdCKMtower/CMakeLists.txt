set(cSMHdCKMtower_SRC
  )

if(ENABLE_TWO_SCALE_SOLVER)
  set(cSMHdCKMtower_SRC
    ${cSMHdCKMtower_SRC}
    cSMHdCKM_cSMHdCKMRHN_two_scale_matching.cpp
    )
endif(ENABLE_TWO_SCALE_SOLVER)

if(ENABLE_SEMI_ANALYTIC_SOLVER)

endif(ENABLE_SEMI_ANALYTIC_SOLVER)

add_library(cSMHdCKMtower ${cSMHdCKMtower_SRC})

if(ENABLE_COMPILER_WARNINGS)
  target_compile_options(cSMHdCKMtower PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${CXX_WARNING_FLAGS}>)
endif(ENABLE_COMPILER_WARNINGS)

target_include_directories(cSMHdCKMtower
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cSMHdCKMtower
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

add_library(FlexibleSUSY::cSMHdCKMtower ALIAS cSMHdCKMtower)

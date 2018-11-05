#include "cSMHdCKM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of Kappa.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_soft_parameters::calc_beta_Kappa_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;

   Eigen::Matrix<std::complex<double>,3,3> beta_Kappa;

   beta_Kappa = (oneOver16PiSqr*(Kappa*(6*traceYdAdjYd + 2*traceYeAdjYe + 6*
                                    traceYuAdjYu - 3*Sqr(g2) + Lambdax)
                               - 1.5*(Ye.transpose()*Ye.conjugate()*Kappa))).
      template cast<std::complex<double> >();

   return beta_Kappa;
}

/**
 * Calculates the 2-loop beta function of Kappa.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_soft_parameters::calc_beta_Kappa_2_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<std::complex<double>,3,3> beta_Kappa;

   beta_Kappa = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();

   return beta_Kappa;
}

/**
 * Calculates the 3-loop beta function of Kappa.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_soft_parameters::calc_beta_Kappa_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Kappa;

   beta_Kappa = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();

   return beta_Kappa;
}

/**
 * Calculates the 4-loop beta function of Kappa.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_soft_parameters::calc_beta_Kappa_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Kappa;

   beta_Kappa = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();

   return beta_Kappa;
}

} // namespace flexiblesusy

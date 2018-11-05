// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Mon 5 Nov 2018 12:48:08

#include "cSMHdCKMRHN_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of Mv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_soft_parameters::calc_beta_Mv_1_loop(const Soft_traces& soft_traces) const
{


   Eigen::Matrix<std::complex<double>,3,3> beta_Mv;

   beta_Mv = (oneOver16PiSqr*(Mv*Yv.conjugate()*Yv.transpose() + Yv*Yv.adjoint(
      )*Mv)).template cast<std::complex<double> >();


   return beta_Mv;
}

/**
 * Calculates the 2-loop beta function of Mv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_soft_parameters::calc_beta_Mv_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<std::complex<double>,3,3> beta_Mv;

   beta_Mv = (twoLoop*(0.075*(-20*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + traceYvAdjYv) + 17*Sqr(g1) + 85*Sqr(g2))*(Mv*Yv.conjugate(
      )*Yv.transpose()) + 0.075*(-20*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + traceYvAdjYv) + 17*Sqr(g1) + 85*Sqr(g2))*(Yv*Yv.adjoint()*
      Mv) - 0.25*(Mv*Yv.conjugate()*Ye.transpose()*Ye.conjugate()*Yv.transpose(
      )) - 0.25*(Mv*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()
      ) - 0.25*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*Mv) + 4*(Yv*Yv.adjoint()*Mv*Yv.
      conjugate()*Yv.transpose()) - 0.25*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*Mv)))
      .template cast<std::complex<double> >();


   return beta_Mv;
}

/**
 * Calculates the 3-loop beta function of Mv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_soft_parameters::calc_beta_Mv_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Mv;

   beta_Mv = ZEROMATRIXCOMPLEX(3,3);


   return beta_Mv;
}

/**
 * Calculates the 4-loop beta function of Mv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_soft_parameters::calc_beta_Mv_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Mv;

   beta_Mv = ZEROMATRIXCOMPLEX(3,3);


   return beta_Mv;
}

} // namespace flexiblesusy

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

// File generated at Mon 5 Nov 2018 12:48:07

#include "cSMHdCKMRHN_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<std::complex<double>,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + traceYvAdjYv - 0.25*Sqr(g1) - 2.25*Sqr(g2) - 8*Sqr(g3)) +
      1.5*(Yd*Yd.adjoint()*Yd) - 1.5*(Yd*Yu.adjoint()*Yu))).template cast<std::
      complex<double> >();


   return beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   Eigen::Matrix<std::complex<double>,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0016666666666666668*Yd*(127*Quad(g1) + 5*Sqr(g1)*(-15
      *(5*traceYdAdjYd + 15*traceYeAdjYe + 17*traceYuAdjYu + 3*traceYvAdjYv) +
      162*Sqr(g2) - 248*Sqr(g3)) + 75*(46*Quad(g2) - 3*Sqr(g2)*(5*(3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + traceYvAdjYv) + 24*Sqr(g3)
      ) + 2*(27*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd + 9*
      traceYeAdjYeYeAdjYe - 2*traceYeAdjYvYvAdjYe + 27*traceYuAdjYuYuAdjYu + 9*
      traceYvAdjYvYvAdjYv + 432*Quad(g3) - 80*(traceYdAdjYd + traceYuAdjYu)*Sqr
      (g3) - 6*Sqr(Lambdax)))) + (-0.75*(9*traceYdAdjYd + 3*traceYeAdjYe + 9*
      traceYuAdjYu + 3*traceYvAdjYv + 8*Lambdax) + 2.3375*Sqr(g1) + 8.4375*Sqr(
      g2) + 16*Sqr(g3))*(Yd*Yd.adjoint()*Yd) + 0.0125*(-79*Sqr(g1) + 45*Sqr(g2)
      + 20*(5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + traceYvAdjYv) -
      64*Sqr(g3)))*(Yd*Yu.adjoint()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd) - Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu - 0.25*(Yd*Yu.adjoint()*Yu*Yd.
      adjoint()*Yd) + 2.75*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).template cast
      <std::complex<double> >();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Yd;

   beta_Yd = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();


   return beta_Yd;
}

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKMRHN_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Yd;

   beta_Yd = ZEROMATRIXCOMPLEX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy

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

// File generated at Mon 5 Nov 2018 12:47:13

#include "cSMHdCKM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<std::complex<double>,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 2.25*Sqr(g1) - 2.25*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()*Ye))).
      template cast<std::complex<double> >();


   return beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<std::complex<double>,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.005*Ye*(1371*Quad(g1) + 5*Sqr(g1)*(25*traceYdAdjYd +
      75*traceYeAdjYe + 85*traceYuAdjYu + 54*Sqr(g2)) - 25*(46*Quad(g2) - 15*(3
      *traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(g2) - 160*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 6*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax)))) + 0.0375*(-20*(9*traceYdAdjYd + 3*traceYeAdjYe + 9*
      traceYuAdjYu + 8*Lambdax) + 129*Sqr(g1) + 225*Sqr(g2))*(Ye*Ye.adjoint()*
      Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).template cast<std::
      complex<double> >();


   return beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Ye;

   beta_Ye = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();


   return beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> cSMHdCKM_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Ye;

   beta_Ye = ZEROMATRIXCOMPLEX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy

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

#include "NMFSU5_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Y10Pr.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> NMFSU5_susy_parameters::calc_beta_Y10Pr_1_loop(const Susy_traces& /* susy_traces */) const
{


   Eigen::Matrix<std::complex<double>,3,3> beta_Y10Pr;

   beta_Y10Pr = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<
      double> >();


   return beta_Y10Pr;
}

/**
 * Calculates the 2-loop beta function of Y10Pr.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> NMFSU5_susy_parameters::calc_beta_Y10Pr_2_loop(const Susy_traces& /* susy_traces */) const
{


   Eigen::Matrix<std::complex<double>,3,3> beta_Y10Pr;

   beta_Y10Pr = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();


   return beta_Y10Pr;
}

/**
 * Calculates the 3-loop beta function of Y10Pr.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> NMFSU5_susy_parameters::calc_beta_Y10Pr_3_loop(const Susy_traces& /* susy_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Y10Pr;

   beta_Y10Pr = (ZEROMATRIXCOMPLEX(3,3)).template cast<std::complex<double> >();


   return beta_Y10Pr;
}

/**
 * Calculates the 4-loop beta function of Y10Pr.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<std::complex<double>,3,3> NMFSU5_susy_parameters::calc_beta_Y10Pr_4_loop(const Susy_traces& /* susy_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<std::complex<double>,3,3> beta_Y10Pr;

   beta_Y10Pr = ZEROMATRIXCOMPLEX(3,3);


   return beta_Y10Pr;
}

} // namespace flexiblesusy

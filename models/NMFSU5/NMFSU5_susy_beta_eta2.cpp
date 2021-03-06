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
 * Calculates the 1-loop beta function of eta2.
 *
 * @return 1-loop beta function
 */
std::complex<double> NMFSU5_susy_parameters::calc_beta_eta2_1_loop(const Susy_traces& /* susy_traces */) const
{


   std::complex<double> beta_eta2;

   beta_eta2 = 0;


   return beta_eta2;
}

/**
 * Calculates the 2-loop beta function of eta2.
 *
 * @return 2-loop beta function
 */
std::complex<double> NMFSU5_susy_parameters::calc_beta_eta2_2_loop(const Susy_traces& /* susy_traces */) const
{


   std::complex<double> beta_eta2;

   beta_eta2 = 0;


   return beta_eta2;
}

/**
 * Calculates the 3-loop beta function of eta2.
 *
 * @return 3-loop beta function
 */
std::complex<double> NMFSU5_susy_parameters::calc_beta_eta2_3_loop(const Susy_traces& /* susy_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   std::complex<double> beta_eta2;

   beta_eta2 = 0;


   return beta_eta2;
}

/**
 * Calculates the 4-loop beta function of eta2.
 *
 * @return 4-loop beta function
 */
std::complex<double> NMFSU5_susy_parameters::calc_beta_eta2_4_loop(const Susy_traces& /* susy_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   std::complex<double> beta_eta2;

   beta_eta2 = 0;


   return beta_eta2;
}

} // namespace flexiblesusy

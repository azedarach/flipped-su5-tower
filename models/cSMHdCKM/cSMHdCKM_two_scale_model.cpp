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

// File generated at Mon 5 Nov 2018 12:47:54

/**
 * @file cSMHdCKM_two_scale_model.cpp
 * @brief implementation of the cSMHdCKM model class
 *
 * Contains the definition of the cSMHdCKM model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated at Mon 5 Nov 2018 12:47:54 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#include "cSMHdCKM_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME cSMHdCKM<Two_scale>

CLASSNAME::cSMHdCKM(const cSMHdCKM_input_parameters& input_)
   : cSMHdCKM_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   cSMHdCKM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   cSMHdCKM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return cSMHdCKM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   cSMHdCKM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   cSMHdCKM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   cSMHdCKM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

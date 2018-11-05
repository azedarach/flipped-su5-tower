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

// File generated at Mon 5 Nov 2018 12:48:09

#include "cSMHdCKMRHN_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd cSMHdCKMRHN_input_parameters::get() const
{
   Eigen::ArrayXd pars(12);

   pars(0) = LambdaIN;
   pars(1) = Qin;
   pars(2) = QEWSB;
   pars(3) = MvInput(0,0);
   pars(4) = MvInput(0,1);
   pars(5) = MvInput(0,2);
   pars(6) = MvInput(1,0);
   pars(7) = MvInput(1,1);
   pars(8) = MvInput(1,2);
   pars(9) = MvInput(2,0);
   pars(10) = MvInput(2,1);
   pars(11) = MvInput(2,2);

   return pars;
}

void cSMHdCKMRHN_input_parameters::set(const Eigen::ArrayXd& pars)
{
   LambdaIN = pars(0);
   Qin = pars(1);
   QEWSB = pars(2);
   MvInput(0,0) = pars(3);
   MvInput(0,1) = pars(4);
   MvInput(0,2) = pars(5);
   MvInput(1,0) = pars(6);
   MvInput(1,1) = pars(7);
   MvInput(1,2) = pars(8);
   MvInput(2,0) = pars(9);
   MvInput(2,1) = pars(10);
   MvInput(2,2) = pars(11);

}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKMRHN_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "QEWSB = " << INPUT(QEWSB) << ", ";
   ostr << "MvInput = " << INPUT(MvInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy

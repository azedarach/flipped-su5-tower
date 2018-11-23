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

#include "cSMHdCKMRHNEFT_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd cSMHdCKMRHNEFT_input_parameters::get() const
{
   Eigen::ArrayXd pars(21);

   pars(0) = LambdaIN;
   pars(1) = Qin;
   pars(2) = sign_delta_mAsq;
   pars(3) = Re(UvInput(0,0));
   pars(4) = Im(UvInput(0,0));
   pars(5) = Re(UvInput(0,1));
   pars(6) = Im(UvInput(0,1));
   pars(7) = Re(UvInput(0,2));
   pars(8) = Im(UvInput(0,2));
   pars(9) = Re(UvInput(1,0));
   pars(10) = Im(UvInput(1,0));
   pars(11) = Re(UvInput(1,1));
   pars(12) = Im(UvInput(1,1));
   pars(13) = Re(UvInput(1,2));
   pars(14) = Im(UvInput(1,2));
   pars(15) = Re(UvInput(2,0));
   pars(16) = Im(UvInput(2,0));
   pars(17) = Re(UvInput(2,1));
   pars(18) = Im(UvInput(2,1));
   pars(19) = Re(UvInput(2,2));
   pars(20) = Im(UvInput(2,2));

   return pars;
}

void cSMHdCKMRHNEFT_input_parameters::set(const Eigen::ArrayXd& pars)
{
   LambdaIN = pars(0);
   Qin = pars(1);
   sign_delta_mAsq = pars(2);
   UvInput(0,0) = std::complex<double>(pars(3), pars(4));
   UvInput(0,1) = std::complex<double>(pars(5), pars(6));
   UvInput(0,2) = std::complex<double>(pars(7), pars(8));
   UvInput(1,0) = std::complex<double>(pars(9), pars(10));
   UvInput(1,1) = std::complex<double>(pars(11), pars(12));
   UvInput(1,2) = std::complex<double>(pars(13), pars(14));
   UvInput(2,0) = std::complex<double>(pars(15), pars(16));
   UvInput(2,1) = std::complex<double>(pars(17), pars(18));
   UvInput(2,2) = std::complex<double>(pars(19), pars(20));
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKMRHNEFT_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "sign_delta_mAsq = " << INPUT(sign_delta_mAsq) << ", ";
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         ostr << "UvInput(" << i << ", " << j << ") = " << INPUT(UvInput(i, j)) << ", ";
      }
   }

   return ostr;
}

} // namespace flexiblesusy

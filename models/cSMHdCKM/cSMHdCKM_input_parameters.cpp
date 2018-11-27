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

// File generated at Mon 5 Nov 2018 12:47:14

#include "cSMHdCKM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd cSMHdCKM_input_parameters::get() const
{
   Eigen::ArrayXd pars(22);

   pars(0) = LambdaIN;
   pars(1) = Qin;
   pars(2) = QEWSB;
   pars(3) = sign_delta_mAsq;
   pars(4) = Re(UvInput(0,0));
   pars(5) = Im(UvInput(0,0));
   pars(6) = Re(UvInput(0,1));
   pars(7) = Im(UvInput(0,1));
   pars(8) = Re(UvInput(0,2));
   pars(9) = Im(UvInput(0,2));
   pars(10) = Re(UvInput(1,0));
   pars(11) = Im(UvInput(1,0));
   pars(12) = Re(UvInput(1,1));
   pars(13) = Im(UvInput(1,1));
   pars(14) = Re(UvInput(1,2));
   pars(15) = Im(UvInput(1,2));
   pars(16) = Re(UvInput(2,0));
   pars(17) = Im(UvInput(2,0));
   pars(18) = Re(UvInput(2,1));
   pars(19) = Im(UvInput(2,1));
   pars(20) = Re(UvInput(2,2));
   pars(21) = Im(UvInput(2,2));

   return pars;
}

void cSMHdCKM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   LambdaIN = pars(0);
   Qin = pars(1);
   QEWSB = pars(2);
   sign_delta_mAsq = pars(3);
   UvInput(0,0) = std::complex<double>(pars(4), pars(5));
   UvInput(0,1) = std::complex<double>(pars(6), pars(7));
   UvInput(0,2) = std::complex<double>(pars(8), pars(9));
   UvInput(1,0) = std::complex<double>(pars(10), pars(11));
   UvInput(1,1) = std::complex<double>(pars(12), pars(13));
   UvInput(1,2) = std::complex<double>(pars(14), pars(15));
   UvInput(2,0) = std::complex<double>(pars(16), pars(17));
   UvInput(2,1) = std::complex<double>(pars(18), pars(19));
   UvInput(2,2) = std::complex<double>(pars(20), pars(21));
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "QEWSB = " << INPUT(QEWSB) << ", ";
   ostr << "sign_delta_mAsq = " << INPUT(sign_delta_mAsq) << ", ";
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         ostr << "UvInput(" << i << ", " << j << ") = " << INPUT(UvInput(i, j)) << ", ";
      }
   }

   return ostr;
}

} // namespace flexiblesusy

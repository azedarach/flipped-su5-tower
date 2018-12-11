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

#include "NMFSU5_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd NMFSU5_input_parameters::get() const
{
   Eigen::ArrayXd pars(20);

   pars(0) = lam1IN;
   pars(1) = lam2IN;
   pars(2) = lam3IN;
   pars(3) = lam3tIN;
   pars(4) = lam4IN;
   pars(5) = lam4tIN;
   pars(6) = lam5IN;
   pars(7) = lam5tIN;
   pars(8) = lam6IN;
   pars(9) = lam6tIN;
   pars(10) = Re(lam7IN);
   pars(11) = Im(lam7IN);
   pars(12) = Re(lam8IN);
   pars(13) = Im(lam8IN);
   pars(14) = Re(eta1IN);
   pars(15) = Im(eta1IN);
   pars(16) = Re(eta2IN);
   pars(17) = Im(eta2IN);
   pars(18) = Re(eta3IN);
   pars(19) = Im(eta3IN);

   return pars;
}

void NMFSU5_input_parameters::set(const Eigen::ArrayXd& pars)
{
   lam1IN = pars(0);
   lam2IN = pars(1);
   lam3IN = pars(2);
   lam3tIN = pars(3);
   lam4IN = pars(4);
   lam4tIN = pars(5);
   lam5IN = pars(6);
   lam5tIN = pars(7);
   lam6IN = pars(8);
   lam6tIN = pars(9);
   lam7IN = std::complex<double>(pars(10), pars(11));
   lam8IN = std::complex<double>(pars(12), pars(13));
   eta1IN = std::complex<double>(pars(14), pars(15));
   eta2IN = std::complex<double>(pars(16), pars(17));
   eta3IN = std::complex<double>(pars(18), pars(19));
}

std::ostream& operator<<(std::ostream& ostr, const NMFSU5_input_parameters& input)
{
   ostr << "lam1IN = " << INPUT(lam1IN) << ", ";
   ostr << "lam2IN = " << INPUT(lam2IN) << ", ";
   ostr << "lam3IN = " << INPUT(lam3IN) << ", ";
   ostr << "lam3tIN = " << INPUT(lam3tIN) << ", ";
   ostr << "lam4IN = " << INPUT(lam4IN) << ", ";
   ostr << "lam4tIN = " << INPUT(lam4tIN) << ", ";
   ostr << "lam5IN = " << INPUT(lam5IN) << ", ";
   ostr << "lam5tIN = " << INPUT(lam5tIN) << ", ";
   ostr << "lam6IN = " << INPUT(lam6IN) << ", ";
   ostr << "lam6tIN = " << INPUT(lam6tIN) << ", ";
   ostr << "lam7IN = " << INPUT(lam7IN) << ", ";
   ostr << "lam8IN = " << INPUT(lam8IN) << ", ";
   ostr << "eta1IN = " << INPUT(eta1IN) << ", ";
   ostr << "eta2IN = " << INPUT(eta2IN) << ", ";
   ostr << "eta3IN = " << INPUT(eta3IN) << ", ";

   return ostr;
}

} // namespace flexiblesusy

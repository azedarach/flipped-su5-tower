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
   Eigen::ArrayXd pars(13);

   pars(0) = LambdaIN;
   pars(1) = Qin;
   pars(2) = QEWSB;
   pars(3) = sign_delta_mAsq;
   pars(4) = UV_theta21;
   pars(5) = UV_theta31;
   pars(6) = UV_theta32;
   pars(7) = UV_phi21;
   pars(8) = UV_phi31;
   pars(9) = UV_phi32;
   pars(10) = UV_chi21;
   pars(11) = UV_chi32;
   pars(12) = UV_gamma;

   return pars;
}

void cSMHdCKM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   LambdaIN = pars(0);
   Qin = pars(1);
   QEWSB = pars(2);
   sign_delta_mAsq = pars(3);
   UV_theta21 = pars(4);
   UV_theta31 = pars(5);
   UV_theta32 = pars(6);
   UV_phi21 = pars(7);
   UV_phi31 = pars(8);
   UV_phi32 = pars(9);
   UV_chi21 = pars(10);
   UV_chi32 = pars(11);
   UV_gamma = pars(12);

}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "QEWSB = " << INPUT(QEWSB) << ", ";
   ostr << "sign_delta_mAsq = " << INPUT(sign_delta_mAsq) << ", ";
   ostr << "UV_theta21 = " << INPUT(UV_theta21) << ", ";
   ostr << "UV_theta31 = " << INPUT(UV_theta31) << ", ";
   ostr << "UV_theta32 = "<< INPUT(UV_theta32) << ", ";
   ostr << "UV_phi21 = " << INPUT(UV_phi21) << ", ";
   ostr << "UV_phi31 = " << INPUT(UV_phi31) << ", ";
   ostr << "UV_phi32 = " << INPUT(UV_phi32) << ", ";
   ostr << "UV_chi21 = " << INPUT(UV_chi21) << ", ";
   ostr << "UV_chi32 = " << INPUT(UV_chi32) << ", ";
   ostr << "UV_gamma = " << INPUT(UV_gamma) << ", ";

   return ostr;
}

} // namespace flexiblesusy

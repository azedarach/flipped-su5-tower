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

#ifndef cSMHdCKMRHNEFT_INPUT_PARAMETERS_H
#define cSMHdCKMRHNEFT_INPUT_PARAMETERS_H

#include "wrappers.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct cSMHdCKMRHNEFT_input_parameters {
   double LambdaIN{};
   double Qin{};
   int sign_delta_mAsq{1};
   double UV_theta21{0.};
   double UV_theta31{0.};
   double UV_theta32{0.};
   double UV_phi21{Pi};
   double UV_phi31{0.};
   double UV_phi32{0.};
   double UV_chi21{0.};
   double UV_chi32{0.};
   double UV_gamma{0.};

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const cSMHdCKMRHNEFT_input_parameters&);

} // namespace flexiblesusy

#endif

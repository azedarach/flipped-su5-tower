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

#ifndef NMFSU5_INPUT_PARAMETERS_H
#define NMFSU5_INPUT_PARAMETERS_H

#include "wrappers.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NMFSU5_input_parameters {
   double lam1IN{};
   double lam2IN{};
   double lam3IN{};
   double lam3tIN{};
   double lam4IN{};
   double lam4tIN{};
   double lam5IN{};
   double lam5tIN{};
   double lam6IN{};
   double lam6tIN{};
   std::complex<double> lam7IN{};
   std::complex<double> lam8IN{};
   std::complex<double> eta1IN{};
   std::complex<double> eta2IN{};
   std::complex<double> eta3IN{};

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const NMFSU5_input_parameters&);

} // namespace flexiblesusy

#endif

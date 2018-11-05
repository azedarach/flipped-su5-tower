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
 * @file cSMHdCKM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 5 Nov 2018 12:47:54 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#ifndef cSMHdCKM_TWO_SCALE_H
#define cSMHdCKM_TWO_SCALE_H

#include "cSMHdCKM_model.hpp"
#include "cSMHdCKM_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class cSMHdCKM<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class cSMHdCKM<Two_scale> : public Model, public cSMHdCKM_mass_eigenstates {
public:
   explicit cSMHdCKM(const cSMHdCKM_input_parameters& input_ = cSMHdCKM_input_parameters());
   cSMHdCKM(const cSMHdCKM&) = default;
   cSMHdCKM(cSMHdCKM&&) = default;
   virtual ~cSMHdCKM() = default;
   cSMHdCKM& operator=(const cSMHdCKM&) = default;
   cSMHdCKM& operator=(cSMHdCKM&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const cSMHdCKM<Two_scale>&);

} // namespace flexiblesusy

#endif

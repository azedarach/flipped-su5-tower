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

// File generated at Mon 5 Nov 2018 12:48:55

/**
 * @file cSMHdCKMRHN_a_muon.hpp
 *
 * This file was generated at Mon 5 Nov 2018 12:48:55 with FlexibleSUSY
 * 2.2.0 and SARAH 4.14.0 .
 */

#ifndef cSMHdCKMRHN_A_MUON_H
#define cSMHdCKMRHN_A_MUON_H

namespace flexiblesusy {
class cSMHdCKMRHN_mass_eigenstates;

namespace cSMHdCKMRHN_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const cSMHdCKMRHN_mass_eigenstates& model);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const cSMHdCKMRHN_mass_eigenstates& model);
} // namespace cSMHdCKMRHN_a_muon
} // namespace flexiblesusy

#endif

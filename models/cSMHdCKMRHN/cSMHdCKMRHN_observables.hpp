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

#ifndef cSMHdCKMRHN_OBSERVABLES_H
#define cSMHdCKMRHN_OBSERVABLES_H

#include <string>
#include <vector>
#include <Eigen/Core>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class cSMHdCKMRHN_mass_eigenstates;
class Physical_input;

struct cSMHdCKMRHN_observables {
   static const int NUMBER_OF_OBSERVABLES = 5;

   cSMHdCKMRHN_observables();
   Eigen::ArrayXd get() const; ///< returns vector of all observables
   static std::vector<std::string> get_names(); ///< returns vector of all observable names
   void clear(); ///< sets all observables to zero
   void set(const Eigen::ArrayXd&); ///< sets all observables from given vector

   double a_muon; ///< a_muon = (g-2)/2 of the muon (calculated with FlexibleSUSY)
   std::complex<double> eff_cp_higgs_photon_photon; ///< effective H-Photon-Photon coupling
   std::complex<double> eff_cp_higgs_gluon_gluon; ///< effective H-Gluon-Gluon coupling

};

cSMHdCKMRHN_observables calculate_observables(
   cSMHdCKMRHN_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&);

cSMHdCKMRHN_observables calculate_observables(
   cSMHdCKMRHN_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&, double scale);

} // namespace flexiblesusy

#endif

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

// File generated at Mon 5 Nov 2018 12:48:53

#ifndef cSMHdCKMRHN_TWO_SCALE_CONVERGENCE_TESTER_H
#define cSMHdCKMRHN_TWO_SCALE_CONVERGENCE_TESTER_H

#include "cSMHdCKMRHN_convergence_tester.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class cSMHdCKMRHN_convergence_tester<Two_scale> : public Convergence_tester_DRbar<cSMHdCKMRHN<Two_scale> > {
public:
   using Scale_getter = Convergence_tester_DRbar<cSMHdCKMRHN<Two_scale>>::Scale_getter;

   cSMHdCKMRHN_convergence_tester(cSMHdCKMRHN<Two_scale>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~cSMHdCKMRHN_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif

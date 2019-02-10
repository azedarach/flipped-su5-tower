#ifndef NMFSU5EFT_MATCHING_H
#define NMFSU5EFT_MATCHING_H

#include "cSMHdCKM_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_mass_eigenstates.hpp"

namespace flexiblesusy {
namespace NMFSU5EFT_matching {

void match_high_to_low_scale_model_tree_level(
   cSMHdCKM_mass_eigenstates&, const cSMHdCKMRHN_mass_eigenstates&);
void match_high_to_low_scale_model(
   cSMHdCKM_mass_eigenstates&, const cSMHdCKMRHN_mass_eigenstates&, int);

void match_low_to_high_scale_model_tree_level(
   cSMHdCKMRHN_mass_eigenstates&, const cSMHdCKM_mass_eigenstates&);
void match_low_to_high_scale_model(
   cSMHdCKMRHN_mass_eigenstates&, const cSMHdCKM_mass_eigenstates&, int);


} // namespace cSMHdCKMRHNEFT_matching
} // namespace flexiblesusy

#endif

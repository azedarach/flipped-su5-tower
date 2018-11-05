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

#include "cSMHdCKM_observables.hpp"
#include "cSMHdCKM_mass_eigenstates.hpp"
#include "cSMHdCKM_a_muon.hpp"
#include "cSMHdCKM_edm.hpp"
#include "cSMHdCKM_effective_couplings.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#ifdef ENABLE_GM2Calc
#include "gm2calc_interface.hpp"
#endif

#define MODEL model
#define AMU a_muon
#define AMUUNCERTAINTY a_muon_uncertainty
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define EFFCPHIGGSPHOTONPHOTON eff_cp_higgs_photon_photon
#define EFFCPHIGGSGLUONGLUON eff_cp_higgs_gluon_gluon
#define EFFCPPSEUDOSCALARPHOTONPHOTON eff_cp_pseudoscalar_photon_photon
#define EFFCPPSEUDOSCALARGLUONGLUON eff_cp_pseudoscalar_gluon_gluon

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

const int cSMHdCKM_observables::NUMBER_OF_OBSERVABLES;

cSMHdCKM_observables::cSMHdCKM_observables()
   : a_muon(0)
   , eff_cp_higgs_photon_photon(0)
   , eff_cp_higgs_gluon_gluon(0)

{
}

Eigen::ArrayXd cSMHdCKM_observables::get() const
{
   Eigen::ArrayXd vec(cSMHdCKM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;
   vec(1) = Re(eff_cp_higgs_photon_photon);
   vec(2) = Im(eff_cp_higgs_photon_photon);
   vec(3) = Re(eff_cp_higgs_gluon_gluon);
   vec(4) = Im(eff_cp_higgs_gluon_gluon);

   return vec;
}

std::vector<std::string> cSMHdCKM_observables::get_names()
{
   std::vector<std::string> names(cSMHdCKM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";
   names[1] = "Re(eff_cp_higgs_photon_photon)";
   names[2] = "Im(eff_cp_higgs_photon_photon)";
   names[3] = "Re(eff_cp_higgs_gluon_gluon)";
   names[4] = "Im(eff_cp_higgs_gluon_gluon)";

   return names;
}

void cSMHdCKM_observables::clear()
{
   a_muon = 0.;
   eff_cp_higgs_photon_photon = std::complex<double>(0.,0.);
   eff_cp_higgs_gluon_gluon = std::complex<double>(0.,0.);

}

void cSMHdCKM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == cSMHdCKM_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);
   eff_cp_higgs_photon_photon = std::complex<double>(vec(1), vec(2));
   eff_cp_higgs_gluon_gluon = std::complex<double>(vec(3), vec(4));

}

cSMHdCKM_observables calculate_observables(cSMHdCKM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const Error& e) {
         model.get_problems().flag_thrown(e.what());
         return cSMHdCKM_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

cSMHdCKM_observables calculate_observables(cSMHdCKM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   cSMHdCKM_observables observables;

   try {
      cSMHdCKM_effective_couplings effective_couplings(model, qedqcd, physical_input);
      effective_couplings.calculate_effective_couplings();

      observables.AMU = cSMHdCKM_a_muon::calculate_a_muon(MODEL);
      observables.EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_CphhVPVP();
      observables.EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_CphhVGVG();
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy

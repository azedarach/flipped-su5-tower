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

// File generated at Mon 5 Nov 2018 12:48:54

#include "cSMHdCKMRHNEFT_two_scale_susy_scale_constraint.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"

#include <cmath>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) input.p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME cSMHdCKMRHN<Two_scale>

cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::cSMHdCKMRHNEFT_susy_scale_constraint(
   cSMHdCKMRHN<Two_scale>* model_, const softsusy::QedQcd& qedqcd_,
   const cSMHdCKMRHNEFT_input_parameters& input_)
   : model(model_)
   , qedqcd(qedqcd_)
   , input(input_)
{
   initialize();
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();




   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   MODEL->solve_ewsb();

}

double cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const cSMHdCKMRHNEFT_input_parameters& cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   return input;
}

cSMHdCKMRHN<Two_scale>* cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<cSMHdCKMRHN<Two_scale>*>(model_);
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
   input = cSMHdCKMRHNEFT_input_parameters();
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   initial_scale_guess = calculate_initial_scale_guess();

   scale = initial_scale_guess;
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MFv = MODELPARAMETER(MFv);

   scale = Cbrt(Abs(MFv(3) * MFv(4) * MFv(5)));


}

double cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::calculate_initial_scale_guess() const
{
   check_model_ptr();

   const auto mu_guess = qedqcd.displayMass(softsusy::mUp);
   const auto mc_guess = qedqcd.displayMass(softsusy::mCharm);
   const auto mt_guess = model->get_thresholds() > 0 && model->get_threshold_corrections().mt > 0 ?
      qedqcd.displayMass(softsusy::mTop) - 30.0 :
      qedqcd.displayPoleMt();
   const auto mv1_guess = qedqcd.displayNeutrinoPoleMass(1);
   const auto mv2_guess = qedqcd.displayNeutrinoPoleMass(2);
   const auto mv3_guess = qedqcd.displayNeutrinoPoleMass(3);
   const auto v = Electroweak_constants::vev;

   Eigen::Matrix<std::complex<double>,3,3> Yu_guess(Eigen::Matrix<std::complex<double>,3,3>::Zero());
   Yu_guess(0,0) = Sqrt(2.) * mu_guess / v;
   Yu_guess(1,1) = Sqrt(2.) * mc_guess / v;
   Yu_guess(2,2) = Sqrt(2.) * mt_guess / v;

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fv_light(Eigen::Matrix<std::complex<double>,3,3>::Zero());
   mass_matrix_Fv_light(0,0) = mv1_guess;
   mass_matrix_Fv_light(1,1) = mv2_guess;
   mass_matrix_Fv_light(2,2) = mv3_guess;

   const Eigen::Matrix<std::complex<double>,3,3> Kappa_guess = 4. * mass_matrix_Fv_light / Sqr(v);

   Eigen::Matrix<std::complex<double>,3,3> Mv_guess(Eigen::Matrix<std::complex<double>,3,3>::Zero());
   calculate_seesaw_Mv(Kappa_guess, Yu_guess.transpose(), Mv_guess);

   return Cbrt(Abs(Mv_guess.determinant()));
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::calculate_seesaw_Mv(
   const Eigen::Matrix<std::complex<double>,3,3>& Kappa, const Eigen::Matrix<std::complex<double>,3,3>& Yv,
   Eigen::Matrix<std::complex<double>,3,3>& Mv) const
{
   Mv = 2. * Yv * Kappa.inverse() * Yv.transpose();
}

void cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy

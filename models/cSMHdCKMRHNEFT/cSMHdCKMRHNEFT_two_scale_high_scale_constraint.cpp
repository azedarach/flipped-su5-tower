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

#include "cSMHdCKMRHNEFT_two_scale_high_scale_constraint.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"
#include "cSMHdCKMRHN_info.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#endif

#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME cSMHdCKMRHN<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
      const auto Vu = MODELPARAMETER(Vu); \
      const auto Ud = MODELPARAMETER(Ud); \
      const auto Uu = MODELPARAMETER(Uu); \
      const auto Ve = MODELPARAMETER(Ve); \
      const auto Ue = MODELPARAMETER(Ue); \
      const auto MAh = MODELPARAMETER(MAh); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(Mu); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(VEVSM1); \
      pars.vu = Re(VEVSM2); \
      pars.mq2 = Re(Vu); \
      pars.md2 = Re(Ud); \
      pars.mu2 = Re(Uu); \
      pars.ml2 = Re(Ve); \
      pars.me2 = Re(Ue); \
      pars.Au(2,2) = Re(TrilinearUp); \
      pars.Ad(2,2) = Re(TrilinearDown); \
      pars.Ae(2,2) = Re(TrilinearLepton); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = 0; \
      pars.M2 = 0; \
      pars.MG = MGluino; \
      pars.MA = MAh; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = 1; \
      const double lambda_3L_uncertainty = 0; \
       \
                                                                        \
      double lambda_3L = 0.;                                            \
                                                                        \
      try {                                                             \
         const bool verbose = false;                                    \
         himalaya::HierarchyCalculator hc(pars, verbose);               \
                                                                        \
         const auto ho = hc.calculateDMh3L(false);                      \
                                                                        \
         lambda_3L =                                                    \
            lambda_3L_eft * (                                           \
               ho.getDLambda(3)                                         \
               + msbar_scheme*ho.getDLambdaDRbarPrimeToMSbarShift(3)    \
               + lambda_3L_uncertainty*ho.getDLambdaUncertainty(3)      \
            );                                                          \
                                                                        \
         VERBOSE_MSG("Himalaya top (hierarchy, Dlambda_3L) = ("         \
                     << ho.getSuitableHierarchy() << ", "               \
                     << lambda_3L <<")");                               \
      } catch (const std::exception& e) {                               \
         model->get_problems().flag_bad_mass(cSMHdCKMRHN_info::hh); \
         WARNING(e.what());                                             \
         VERBOSE_MSG(pars);                                             \
      }                                                                 \
                                                                        \
      return lambda_3L;                                                 \
   }()
#else
#define FSHimalayaMh23L [] () {                                         \
      throw HimalayaError("The 3-loop corrections to lambda "           \
                          "require Himalaya 2.0.0 (or higher), but "    \
                          "FlexibleSUSY has not been configured with "  \
                          "this Himalaya version!");                    \
      return 0.;                                                        \
   }()
#endif

cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::cSMHdCKMRHNEFT_high_scale_constraint(
   cSMHdCKMRHN<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();




   update_scale();

   const auto LambdaIN = INPUTPARAMETER(LambdaIN);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_Lambdax(Re(LambdaIN));
   MODEL->set_Yv((Yu.transpose()).template cast<std::complex<double> >());

   Eigen::Matrix<std::complex<double>,3,3> Ud;
   Eigen::Matrix<std::complex<double>,3,3> Vd;
   Eigen::Array<double,3,1> Yd_diag;
   fs_svd(MODELPARAMETER(Yd), Yd_diag, Ud, Vd);

   MODEL->set_Yd(Vd.transpose() * Yd_diag.matrix().asDiagonal() * Vd);

   check_non_perturbative();
}

bool cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Yu = MODELPARAMETER(Yu);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g3);
   }
   if (MaxAbsValue(Lambdax) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::Lambdax, MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::Lambdax);
   }
   if (MaxAbsValue(Re(Yd(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_0, MaxAbsValue(Re(Yd(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_0);
   }

   if (MaxAbsValue(Im(Yd(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_0, MaxAbsValue(Im(Yd(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_0);
   }

   if (MaxAbsValue(Re(Yd(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_1, MaxAbsValue(Re(Yd(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_1);
   }

   if (MaxAbsValue(Im(Yd(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_1, MaxAbsValue(Im(Yd(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_1);
   }

   if (MaxAbsValue(Re(Yd(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_2, MaxAbsValue(Re(Yd(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_2);
   }

   if (MaxAbsValue(Im(Yd(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_2, MaxAbsValue(Im(Yd(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_2);
   }

   if (MaxAbsValue(Re(Yd(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_0, MaxAbsValue(Re(Yd(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_0);
   }

   if (MaxAbsValue(Im(Yd(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_0, MaxAbsValue(Im(Yd(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_0);
   }

   if (MaxAbsValue(Re(Yd(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_1, MaxAbsValue(Re(Yd(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_1);
   }

   if (MaxAbsValue(Im(Yd(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_1, MaxAbsValue(Im(Yd(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_1);
   }

   if (MaxAbsValue(Re(Yd(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_2, MaxAbsValue(Re(Yd(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_2);
   }

   if (MaxAbsValue(Im(Yd(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_2, MaxAbsValue(Im(Yd(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_2);
   }

   if (MaxAbsValue(Re(Yd(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_0, MaxAbsValue(Re(Yd(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_0);
   }

   if (MaxAbsValue(Im(Yd(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_0, MaxAbsValue(Im(Yd(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_0);
   }

   if (MaxAbsValue(Re(Yd(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_1, MaxAbsValue(Re(Yd(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_1);
   }

   if (MaxAbsValue(Im(Yd(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_1, MaxAbsValue(Im(Yd(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_1);
   }

   if (MaxAbsValue(Re(Yd(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_2, MaxAbsValue(Re(Yd(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_2);
   }

   if (MaxAbsValue(Im(Yd(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_2, MaxAbsValue(Im(Yd(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_2);
   }
   if (MaxAbsValue(Re(Ye(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_0, MaxAbsValue(Re(Ye(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_0);
   }

   if (MaxAbsValue(Im(Ye(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_0, MaxAbsValue(Im(Ye(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_0);
   }

   if (MaxAbsValue(Re(Ye(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_1, MaxAbsValue(Re(Ye(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_1);
   }

   if (MaxAbsValue(Im(Ye(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_1, MaxAbsValue(Im(Ye(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_1);
   }

   if (MaxAbsValue(Re(Ye(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_2, MaxAbsValue(Re(Ye(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_2);
   }

   if (MaxAbsValue(Im(Ye(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_2, MaxAbsValue(Im(Ye(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_2);
   }

   if (MaxAbsValue(Re(Ye(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_0, MaxAbsValue(Re(Ye(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_0);
   }

   if (MaxAbsValue(Im(Ye(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_0, MaxAbsValue(Im(Ye(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_0);
   }

   if (MaxAbsValue(Re(Ye(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_1, MaxAbsValue(Re(Ye(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_1);
   }

   if (MaxAbsValue(Im(Ye(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_1, MaxAbsValue(Im(Ye(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_1);
   }

   if (MaxAbsValue(Re(Ye(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_2, MaxAbsValue(Re(Ye(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_2);
   }

   if (MaxAbsValue(Im(Ye(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_2, MaxAbsValue(Im(Ye(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_2);
   }

   if (MaxAbsValue(Re(Ye(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_0, MaxAbsValue(Re(Ye(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_0);
   }

   if (MaxAbsValue(Im(Ye(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_0, MaxAbsValue(Im(Ye(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_0);
   }

   if (MaxAbsValue(Re(Ye(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_1, MaxAbsValue(Re(Ye(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_1);
   }

   if (MaxAbsValue(Im(Ye(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_1, MaxAbsValue(Im(Ye(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_1);
   }

   if (MaxAbsValue(Re(Ye(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_2, MaxAbsValue(Re(Ye(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_2);
   }

   if (MaxAbsValue(Im(Ye(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_2, MaxAbsValue(Im(Ye(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_2);
   }
   if (MaxAbsValue(Re(Yv(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_0, MaxAbsValue(Re(Yv(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_0);
   }

   if (MaxAbsValue(Im(Yv(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_0, MaxAbsValue(Im(Yv(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_0);
   }

   if (MaxAbsValue(Re(Yv(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_1, MaxAbsValue(Re(Yv(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_1);
   }

   if (MaxAbsValue(Im(Yv(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_1, MaxAbsValue(Im(Yv(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_1);
   }

   if (MaxAbsValue(Re(Yv(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_2, MaxAbsValue(Re(Yv(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_2);
   }

   if (MaxAbsValue(Im(Yv(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_2, MaxAbsValue(Im(Yv(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_2);
   }

   if (MaxAbsValue(Re(Yv(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_0, MaxAbsValue(Re(Yv(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_0);
   }

   if (MaxAbsValue(Im(Yv(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_0, MaxAbsValue(Im(Yv(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_0);
   }

   if (MaxAbsValue(Re(Yv(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_1, MaxAbsValue(Re(Yv(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_1);
   }

   if (MaxAbsValue(Im(Yv(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_1, MaxAbsValue(Im(Yv(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_1);
   }

   if (MaxAbsValue(Re(Yv(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_2, MaxAbsValue(Re(Yv(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_2);
   }

   if (MaxAbsValue(Im(Yv(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_2, MaxAbsValue(Im(Yv(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_2);
   }

   if (MaxAbsValue(Re(Yv(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_0, MaxAbsValue(Re(Yv(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_0);
   }

   if (MaxAbsValue(Im(Yv(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_0, MaxAbsValue(Im(Yv(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_0);
   }

   if (MaxAbsValue(Re(Yv(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_1, MaxAbsValue(Re(Yv(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_1);
   }

   if (MaxAbsValue(Im(Yv(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_1, MaxAbsValue(Im(Yv(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_1);
   }

   if (MaxAbsValue(Re(Yv(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_2, MaxAbsValue(Re(Yv(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_2);
   }

   if (MaxAbsValue(Im(Yv(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_2, MaxAbsValue(Im(Yv(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_2);
   }
   if (MaxAbsValue(Re(Yu(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_0, MaxAbsValue(Re(Yu(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_0);
   }

   if (MaxAbsValue(Im(Yu(0,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_0, MaxAbsValue(Im(Yu(0,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_0);
   }

   if (MaxAbsValue(Re(Yu(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_1, MaxAbsValue(Re(Yu(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_1);
   }

   if (MaxAbsValue(Im(Yu(0,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_1, MaxAbsValue(Im(Yu(0,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_1);
   }

   if (MaxAbsValue(Re(Yu(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_2, MaxAbsValue(Re(Yu(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_2);
   }

   if (MaxAbsValue(Im(Yu(0,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_2, MaxAbsValue(Im(Yu(0,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_2);
   }

   if (MaxAbsValue(Re(Yu(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_0, MaxAbsValue(Re(Yu(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_0);
   }

   if (MaxAbsValue(Im(Yu(1,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_0, MaxAbsValue(Im(Yu(1,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_0);
   }

   if (MaxAbsValue(Re(Yu(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_1, MaxAbsValue(Re(Yu(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_1);
   }

   if (MaxAbsValue(Im(Yu(1,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_1, MaxAbsValue(Im(Yu(1,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_1);
   }

   if (MaxAbsValue(Re(Yu(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_2, MaxAbsValue(Re(Yu(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_2);
   }

   if (MaxAbsValue(Im(Yu(1,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_2, MaxAbsValue(Im(Yu(1,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_2);
   }

   if (MaxAbsValue(Re(Yu(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_0, MaxAbsValue(Re(Yu(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_0);
   }

   if (MaxAbsValue(Im(Yu(2,0))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_0, MaxAbsValue(Im(Yu(2,0))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_0);
   }

   if (MaxAbsValue(Re(Yu(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_1, MaxAbsValue(Re(Yu(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_1);
   }

   if (MaxAbsValue(Im(Yu(2,1))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_1, MaxAbsValue(Im(Yu(2,1))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_1);
   }

   if (MaxAbsValue(Re(Yu(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_2, MaxAbsValue(Re(Yu(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_2);
   }

   if (MaxAbsValue(Im(Yu(2,2))) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_2, MaxAbsValue(Im(Yu(2,2))), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_2);
   }


   return problem;
}

double cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const cSMHdCKMRHN_input_parameters& cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

cSMHdCKMRHN<Two_scale>* cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<cSMHdCKMRHN<Two_scale>*>(model_);
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


}

void cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy

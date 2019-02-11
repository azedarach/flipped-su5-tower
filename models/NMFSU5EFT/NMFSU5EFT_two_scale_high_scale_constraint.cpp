#include "NMFSU5EFT_two_scale_high_scale_constraint.hpp"
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
#include "wittens_loop.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#endif

#include <Eigen/Dense>

#include <cmath>
#include <cerrno>
#include <cstring>

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

NMFSU5EFT_high_scale_constraint<Two_scale>::NMFSU5EFT_high_scale_constraint(
   cSMHdCKMRHN<Two_scale>* model_, const NMFSU5EFT_input_parameters& input_)
   : model(model_)
   , input(input_)
{
   initialize();
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   const auto Lambda1IN = INPUTPARAMETER(Lambda1IN);
   const auto Lambda2IN = INPUTPARAMETER(Lambda2IN);
   const auto Lambda3IN = INPUTPARAMETER(Lambda3IN);
   const auto Lambda3tIN = INPUTPARAMETER(Lambda3tIN);
   const auto Lambda4IN = INPUTPARAMETER(Lambda4IN);
   const auto Lambda4tIN = INPUTPARAMETER(Lambda4tIN);
   const auto Lambda5IN = INPUTPARAMETER(Lambda5IN);
   const auto Lambda5tIN = INPUTPARAMETER(Lambda5tIN);
   const auto Lambda6IN = INPUTPARAMETER(Lambda6IN);
   const auto Lambda6tIN = INPUTPARAMETER(Lambda6tIN);
   const auto Lambda7IN = INPUTPARAMETER(Lambda7IN);
   const auto Lambda8IN = INPUTPARAMETER(Lambda8IN);
   const auto Eta1IN = INPUTPARAMETER(Eta1IN);
   const auto Eta2IN = INPUTPARAMETER(Eta2IN);
   const auto Eta3IN = INPUTPARAMETER(Eta3IN);

   HIGHSCALEMODEL->set_lam1(Lambda1IN);
   HIGHSCALEMODEL->set_lam2(Lambda2IN);
   HIGHSCALEMODEL->set_lam3(Lambda3IN);
   HIGHSCALEMODEL->set_lam3t(Lambda3tIN);
   HIGHSCALEMODEL->set_lam4(Lambda4IN);
   HIGHSCALEMODEL->set_lam4t(Lambda4tIN);
   HIGHSCALEMODEL->set_lam5(Lambda5IN);
   HIGHSCALEMODEL->set_lam5t(Lambda5tIN);
   HIGHSCALEMODEL->set_lam6(Lambda6IN);
   HIGHSCALEMODEL->set_lam6t(Lambda6tIN);
   HIGHSCALEMODEL->set_lam7(Lambda7IN);
   HIGHSCALEMODEL->set_lam8(Lambda8IN);
   HIGHSCALEMODEL->set_eta1(Eta1IN);
   HIGHSCALEMODEL->set_eta2(Eta2IN);
   HIGHSCALEMODEL->set_eta3(Eta3IN);

   HIGHSCALEMODEL->calculate_DRbar_masses();

   calculate_threshold_corrections();

   calculate_DRbar_yukawa_couplings();

   update_scale();

   check_non_perturbative();
   check_high_scale_non_perturbative();
}

bool NMFSU5EFT_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   const double perturbativity_bound = 3.5449077018110318;
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yv = MODELPARAMETER(Yv);
   const auto Yu = MODELPARAMETER(Yu);

   if (MaxAbsValue(g1) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g1, MaxAbsValue(g1), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g1);
   }
   if (MaxAbsValue(g2) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g2, MaxAbsValue(g2), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g2);
   }
   if (MaxAbsValue(g3) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::g3, MaxAbsValue(g3), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::g3);
   }
   if (MaxAbsValue(Lambdax) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::Lambdax, MaxAbsValue(Lambdax), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::Lambdax);
   }
   if (MaxAbsValue(Re(Yd(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_0, MaxAbsValue(Re(Yd(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_0);
   }

   if (MaxAbsValue(Im(Yd(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_0, MaxAbsValue(Im(Yd(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_0);
   }

   if (MaxAbsValue(Re(Yd(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_1, MaxAbsValue(Re(Yd(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_1);
   }

   if (MaxAbsValue(Im(Yd(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_1, MaxAbsValue(Im(Yd(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_1);
   }

   if (MaxAbsValue(Re(Yd(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_2, MaxAbsValue(Re(Yd(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd0_2);
   }

   if (MaxAbsValue(Im(Yd(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_2, MaxAbsValue(Im(Yd(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd0_2);
   }

   if (MaxAbsValue(Re(Yd(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_0, MaxAbsValue(Re(Yd(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_0);
   }

   if (MaxAbsValue(Im(Yd(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_0, MaxAbsValue(Im(Yd(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_0);
   }

   if (MaxAbsValue(Re(Yd(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_1, MaxAbsValue(Re(Yd(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_1);
   }

   if (MaxAbsValue(Im(Yd(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_1, MaxAbsValue(Im(Yd(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_1);
   }

   if (MaxAbsValue(Re(Yd(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_2, MaxAbsValue(Re(Yd(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd1_2);
   }

   if (MaxAbsValue(Im(Yd(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_2, MaxAbsValue(Im(Yd(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd1_2);
   }

   if (MaxAbsValue(Re(Yd(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_0, MaxAbsValue(Re(Yd(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_0);
   }

   if (MaxAbsValue(Im(Yd(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_0, MaxAbsValue(Im(Yd(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_0);
   }

   if (MaxAbsValue(Re(Yd(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_1, MaxAbsValue(Re(Yd(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_1);
   }

   if (MaxAbsValue(Im(Yd(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_1, MaxAbsValue(Im(Yd(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_1);
   }

   if (MaxAbsValue(Re(Yd(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_2, MaxAbsValue(Re(Yd(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYd2_2);
   }

   if (MaxAbsValue(Im(Yd(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_2, MaxAbsValue(Im(Yd(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYd2_2);
   }
   if (MaxAbsValue(Re(Ye(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_0, MaxAbsValue(Re(Ye(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_0);
   }

   if (MaxAbsValue(Im(Ye(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_0, MaxAbsValue(Im(Ye(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_0);
   }

   if (MaxAbsValue(Re(Ye(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_1, MaxAbsValue(Re(Ye(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_1);
   }

   if (MaxAbsValue(Im(Ye(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_1, MaxAbsValue(Im(Ye(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_1);
   }

   if (MaxAbsValue(Re(Ye(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_2, MaxAbsValue(Re(Ye(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe0_2);
   }

   if (MaxAbsValue(Im(Ye(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_2, MaxAbsValue(Im(Ye(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe0_2);
   }

   if (MaxAbsValue(Re(Ye(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_0, MaxAbsValue(Re(Ye(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_0);
   }

   if (MaxAbsValue(Im(Ye(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_0, MaxAbsValue(Im(Ye(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_0);
   }

   if (MaxAbsValue(Re(Ye(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_1, MaxAbsValue(Re(Ye(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_1);
   }

   if (MaxAbsValue(Im(Ye(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_1, MaxAbsValue(Im(Ye(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_1);
   }

   if (MaxAbsValue(Re(Ye(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_2, MaxAbsValue(Re(Ye(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe1_2);
   }

   if (MaxAbsValue(Im(Ye(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_2, MaxAbsValue(Im(Ye(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe1_2);
   }

   if (MaxAbsValue(Re(Ye(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_0, MaxAbsValue(Re(Ye(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_0);
   }

   if (MaxAbsValue(Im(Ye(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_0, MaxAbsValue(Im(Ye(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_0);
   }

   if (MaxAbsValue(Re(Ye(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_1, MaxAbsValue(Re(Ye(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_1);
   }

   if (MaxAbsValue(Im(Ye(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_1, MaxAbsValue(Im(Ye(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_1);
   }

   if (MaxAbsValue(Re(Ye(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_2, MaxAbsValue(Re(Ye(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYe2_2);
   }

   if (MaxAbsValue(Im(Ye(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_2, MaxAbsValue(Im(Ye(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYe2_2);
   }
   if (MaxAbsValue(Re(Yv(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_0, MaxAbsValue(Re(Yv(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_0);
   }

   if (MaxAbsValue(Im(Yv(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_0, MaxAbsValue(Im(Yv(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_0);
   }

   if (MaxAbsValue(Re(Yv(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_1, MaxAbsValue(Re(Yv(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_1);
   }

   if (MaxAbsValue(Im(Yv(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_1, MaxAbsValue(Im(Yv(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_1);
   }

   if (MaxAbsValue(Re(Yv(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_2, MaxAbsValue(Re(Yv(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv0_2);
   }

   if (MaxAbsValue(Im(Yv(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_2, MaxAbsValue(Im(Yv(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv0_2);
   }

   if (MaxAbsValue(Re(Yv(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_0, MaxAbsValue(Re(Yv(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_0);
   }

   if (MaxAbsValue(Im(Yv(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_0, MaxAbsValue(Im(Yv(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_0);
   }

   if (MaxAbsValue(Re(Yv(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_1, MaxAbsValue(Re(Yv(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_1);
   }

   if (MaxAbsValue(Im(Yv(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_1, MaxAbsValue(Im(Yv(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_1);
   }

   if (MaxAbsValue(Re(Yv(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_2, MaxAbsValue(Re(Yv(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv1_2);
   }

   if (MaxAbsValue(Im(Yv(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_2, MaxAbsValue(Im(Yv(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv1_2);
   }

   if (MaxAbsValue(Re(Yv(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_0, MaxAbsValue(Re(Yv(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_0);
   }

   if (MaxAbsValue(Im(Yv(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_0, MaxAbsValue(Im(Yv(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_0);
   }

   if (MaxAbsValue(Re(Yv(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_1, MaxAbsValue(Re(Yv(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_1);
   }

   if (MaxAbsValue(Im(Yv(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_1, MaxAbsValue(Im(Yv(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_1);
   }

   if (MaxAbsValue(Re(Yv(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_2, MaxAbsValue(Re(Yv(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYv2_2);
   }

   if (MaxAbsValue(Im(Yv(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_2, MaxAbsValue(Im(Yv(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYv2_2);
   }
   if (MaxAbsValue(Re(Yu(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_0, MaxAbsValue(Re(Yu(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_0);
   }

   if (MaxAbsValue(Im(Yu(0,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_0, MaxAbsValue(Im(Yu(0,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_0);
   }

   if (MaxAbsValue(Re(Yu(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_1, MaxAbsValue(Re(Yu(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_1);
   }

   if (MaxAbsValue(Im(Yu(0,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_1, MaxAbsValue(Im(Yu(0,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_1);
   }

   if (MaxAbsValue(Re(Yu(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_2, MaxAbsValue(Re(Yu(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu0_2);
   }

   if (MaxAbsValue(Im(Yu(0,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_2, MaxAbsValue(Im(Yu(0,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu0_2);
   }

   if (MaxAbsValue(Re(Yu(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_0, MaxAbsValue(Re(Yu(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_0);
   }

   if (MaxAbsValue(Im(Yu(1,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_0, MaxAbsValue(Im(Yu(1,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_0);
   }

   if (MaxAbsValue(Re(Yu(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_1, MaxAbsValue(Re(Yu(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_1);
   }

   if (MaxAbsValue(Im(Yu(1,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_1, MaxAbsValue(Im(Yu(1,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_1);
   }

   if (MaxAbsValue(Re(Yu(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_2, MaxAbsValue(Re(Yu(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu1_2);
   }

   if (MaxAbsValue(Im(Yu(1,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_2, MaxAbsValue(Im(Yu(1,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu1_2);
   }

   if (MaxAbsValue(Re(Yu(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_0, MaxAbsValue(Re(Yu(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_0);
   }

   if (MaxAbsValue(Im(Yu(2,0))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_0, MaxAbsValue(Im(Yu(2,0))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_0);
   }

   if (MaxAbsValue(Re(Yu(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_1, MaxAbsValue(Re(Yu(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_1);
   }

   if (MaxAbsValue(Im(Yu(2,1))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_1, MaxAbsValue(Im(Yu(2,1))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_1);
   }

   if (MaxAbsValue(Re(Yu(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_2, MaxAbsValue(Re(Yu(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ReYu2_2);
   }

   if (MaxAbsValue(Im(Yu(2,2))) > perturbativity_bound) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_2, MaxAbsValue(Im(Yu(2,2))), model->get_scale(), perturbativity_bound);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(cSMHdCKMRHN_info::ImYu2_2);
   }


   return problem;
}

bool NMFSU5EFT_high_scale_constraint<Two_scale>::check_high_scale_non_perturbative()
{
   const double perturbativity_bound = 3.5449077018110318;
   bool problem = false;

   const auto g5 = HIGHSCALEMODELPARAMETER(g5);
   const auto gX = HIGHSCALEMODELPARAMETER(gX);
   const auto lam1 = HIGHSCALEMODELPARAMETER(lam1);
   const auto lam2 = HIGHSCALEMODELPARAMETER(lam2);
   const auto lam3 = HIGHSCALEMODELPARAMETER(lam3);
   const auto lam3t = HIGHSCALEMODELPARAMETER(lam3t);
   const auto lam4 = HIGHSCALEMODELPARAMETER(lam4);
   const auto lam4t = HIGHSCALEMODELPARAMETER(lam4t);
   const auto lam5 = HIGHSCALEMODELPARAMETER(lam5);
   const auto lam5t = HIGHSCALEMODELPARAMETER(lam5t);
   const auto lam6 = HIGHSCALEMODELPARAMETER(lam6);
   const auto lam6t = HIGHSCALEMODELPARAMETER(lam6t);
   const auto lam7 = HIGHSCALEMODELPARAMETER(lam7);
   const auto lam8 = HIGHSCALEMODELPARAMETER(lam8);
   const auto eta1 = HIGHSCALEMODELPARAMETER(eta1);
   const auto eta2 = HIGHSCALEMODELPARAMETER(eta2);
   const auto eta3 = HIGHSCALEMODELPARAMETER(eta3);
   const auto Y10 = HIGHSCALEMODELPARAMETER(Y10);
   const auto Y10Pr = HIGHSCALEMODELPARAMETER(Y10Pr);
   const auto Y5 = HIGHSCALEMODELPARAMETER(Y5);
   const auto Y5Pr = HIGHSCALEMODELPARAMETER(Y5Pr);
   const auto Y1 = HIGHSCALEMODELPARAMETER(Y1);
   const auto Y1Pr = HIGHSCALEMODELPARAMETER(Y1Pr);

   if (MaxAbsValue(g5) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::g5, MaxAbsValue(g5),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::g5);
   }

   if (MaxAbsValue(gX) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::gX, MaxAbsValue(gX),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::gX);
   }

   if (MaxAbsValue(lam1) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam1, MaxAbsValue(lam1),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam1);
   }

   if (MaxAbsValue(lam2) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam2, MaxAbsValue(lam2),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam2);
   }

   if (MaxAbsValue(lam3) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam3, MaxAbsValue(lam3),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam3);
   }

   if (MaxAbsValue(lam3t) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam3t, MaxAbsValue(lam3t),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam3t);
   }

   if (MaxAbsValue(lam4) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam4, MaxAbsValue(lam4),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam4);
   }

   if (MaxAbsValue(lam4t) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam4t, MaxAbsValue(lam4t),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam4t);
   }

   if (MaxAbsValue(lam5) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam5, MaxAbsValue(lam5),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam5);
   }

   if (MaxAbsValue(lam5t) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam5t, MaxAbsValue(lam5t),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam5t);
   }

   if (MaxAbsValue(lam6) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam6, MaxAbsValue(lam6),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam6);
   }

   if (MaxAbsValue(lam6t) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::lam6t, MaxAbsValue(lam6t),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::lam6t);
   }

   if (MaxAbsValue(Re(lam7)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Relam7, MaxAbsValue(Re(lam7)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Relam7);
   }

   if (MaxAbsValue(Im(lam7)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Imlam7, MaxAbsValue(Im(lam7)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Imlam7);
   }

   if (MaxAbsValue(Re(lam8)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Relam8, MaxAbsValue(Re(lam8)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Relam8);
   }

   if (MaxAbsValue(Im(lam8)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Imlam8, MaxAbsValue(Im(lam8)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Imlam8);
   }

   if (MaxAbsValue(Re(eta1)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Reeta1, MaxAbsValue(Re(eta1)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Reeta1);
   }

   if (MaxAbsValue(Im(eta1)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Imeta1, MaxAbsValue(Im(eta1)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Imeta1);
   }

   if (MaxAbsValue(Re(eta2)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Reeta2, MaxAbsValue(Re(eta2)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Reeta2);
   }

   if (MaxAbsValue(Im(eta2)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Imeta2, MaxAbsValue(Im(eta2)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Imeta2);
   }

   if (MaxAbsValue(Re(eta3)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Reeta3, MaxAbsValue(Re(eta3)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Reeta3);
   }

   if (MaxAbsValue(Im(eta3)) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::Imeta3, MaxAbsValue(Im(eta3)),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::Imeta3);
   }

   if (MaxAbsValue(Re(Y10(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY100_0, MaxAbsValue(Re(Y10(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY100_0);
   }

   if (MaxAbsValue(Im(Y10(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY100_0, MaxAbsValue(Im(Y10(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY100_0);
   }

   if (MaxAbsValue(Re(Y10(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY100_1, MaxAbsValue(Re(Y10(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY100_1);
   }

   if (MaxAbsValue(Im(Y10(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY100_1, MaxAbsValue(Im(Y10(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY100_1);
   }

   if (MaxAbsValue(Re(Y10(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY100_2, MaxAbsValue(Re(Y10(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY100_2);
   }

   if (MaxAbsValue(Im(Y10(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY100_2, MaxAbsValue(Im(Y10(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY100_2);
   }

   if (MaxAbsValue(Re(Y10(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY101_0, MaxAbsValue(Re(Y10(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY101_0);
   }

   if (MaxAbsValue(Im(Y10(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY101_0, MaxAbsValue(Im(Y10(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY101_0);
   }

   if (MaxAbsValue(Re(Y10(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY101_1, MaxAbsValue(Re(Y10(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY101_1);
   }

   if (MaxAbsValue(Im(Y10(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY101_1, MaxAbsValue(Im(Y10(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY101_1);
   }

   if (MaxAbsValue(Re(Y10(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY101_2, MaxAbsValue(Re(Y10(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY101_2);
   }

   if (MaxAbsValue(Im(Y10(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY101_2, MaxAbsValue(Im(Y10(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY101_2);
   }

   if (MaxAbsValue(Re(Y10(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY102_0, MaxAbsValue(Re(Y10(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY102_0);
   }

   if (MaxAbsValue(Im(Y10(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY102_0, MaxAbsValue(Im(Y10(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY102_0);
   }

   if (MaxAbsValue(Re(Y10(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY102_1, MaxAbsValue(Re(Y10(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY102_1);
   }

   if (MaxAbsValue(Im(Y10(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY102_1, MaxAbsValue(Im(Y10(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY102_1);
   }

   if (MaxAbsValue(Re(Y10(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY102_2, MaxAbsValue(Re(Y10(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY102_2);
   }

   if (MaxAbsValue(Im(Y10(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY102_2, MaxAbsValue(Im(Y10(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY102_2);
   }

   if (MaxAbsValue(Re(Y10Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_0, MaxAbsValue(Re(Y10Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_0);
   }

   if (MaxAbsValue(Im(Y10Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_0, MaxAbsValue(Im(Y10Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_0);
   }

   if (MaxAbsValue(Re(Y10Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_1, MaxAbsValue(Re(Y10Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_1);
   }

   if (MaxAbsValue(Im(Y10Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_1, MaxAbsValue(Im(Y10Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_1);
   }

   if (MaxAbsValue(Re(Y10Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_2, MaxAbsValue(Re(Y10Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr0_2);
   }

   if (MaxAbsValue(Im(Y10Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_2, MaxAbsValue(Im(Y10Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr0_2);
   }

   if (MaxAbsValue(Re(Y10Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_0, MaxAbsValue(Re(Y10Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_0);
   }

   if (MaxAbsValue(Im(Y10Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_0, MaxAbsValue(Im(Y10Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_0);
   }

   if (MaxAbsValue(Re(Y10Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_1, MaxAbsValue(Re(Y10Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_1);
   }

   if (MaxAbsValue(Im(Y10Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_1, MaxAbsValue(Im(Y10Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_1);
   }

   if (MaxAbsValue(Re(Y10Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_2, MaxAbsValue(Re(Y10Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr1_2);
   }

   if (MaxAbsValue(Im(Y10Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_2, MaxAbsValue(Im(Y10Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr1_2);
   }

   if (MaxAbsValue(Re(Y10Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_0, MaxAbsValue(Re(Y10Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_0);
   }

   if (MaxAbsValue(Im(Y10Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_0, MaxAbsValue(Im(Y10Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_0);
   }

   if (MaxAbsValue(Re(Y10Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_1, MaxAbsValue(Re(Y10Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_1);
   }

   if (MaxAbsValue(Im(Y10Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_1, MaxAbsValue(Im(Y10Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_1);
   }

   if (MaxAbsValue(Re(Y10Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_2, MaxAbsValue(Re(Y10Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10Pr2_2);
   }

   if (MaxAbsValue(Im(Y10Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_2, MaxAbsValue(Im(Y10Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10Pr2_2);
   }

   if (MaxAbsValue(Re(Y5(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY50_0, MaxAbsValue(Re(Y5(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY50_0);
   }

   if (MaxAbsValue(Im(Y5(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY50_0, MaxAbsValue(Im(Y5(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY50_0);
   }

   if (MaxAbsValue(Re(Y5(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY50_1, MaxAbsValue(Re(Y5(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY50_1);
   }

   if (MaxAbsValue(Im(Y5(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY50_1, MaxAbsValue(Im(Y5(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY50_1);
   }

   if (MaxAbsValue(Re(Y5(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY50_2, MaxAbsValue(Re(Y5(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY50_2);
   }

   if (MaxAbsValue(Im(Y5(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY50_2, MaxAbsValue(Im(Y5(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY50_2);
   }

   if (MaxAbsValue(Re(Y5(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY51_0, MaxAbsValue(Re(Y5(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY51_0);
   }

   if (MaxAbsValue(Im(Y5(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY51_0, MaxAbsValue(Im(Y5(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY51_0);
   }

   if (MaxAbsValue(Re(Y5(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY51_1, MaxAbsValue(Re(Y5(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY51_1);
   }

   if (MaxAbsValue(Im(Y5(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY51_1, MaxAbsValue(Im(Y5(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY51_1);
   }

   if (MaxAbsValue(Re(Y5(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY51_2, MaxAbsValue(Re(Y5(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY51_2);
   }

   if (MaxAbsValue(Im(Y5(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY51_2, MaxAbsValue(Im(Y5(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY51_2);
   }

   if (MaxAbsValue(Re(Y5(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY52_0, MaxAbsValue(Re(Y5(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY52_0);
   }

   if (MaxAbsValue(Im(Y5(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY52_0, MaxAbsValue(Im(Y5(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY52_0);
   }

   if (MaxAbsValue(Re(Y5(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY52_1, MaxAbsValue(Re(Y5(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY52_1);
   }

   if (MaxAbsValue(Im(Y5(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY52_1, MaxAbsValue(Im(Y5(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY52_1);
   }

   if (MaxAbsValue(Re(Y5(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY52_2, MaxAbsValue(Re(Y5(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY52_2);
   }

   if (MaxAbsValue(Im(Y5(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY52_2, MaxAbsValue(Im(Y5(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY52_2);
   }

   if (MaxAbsValue(Re(Y5Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_0, MaxAbsValue(Re(Y5Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_0);
   }

   if (MaxAbsValue(Im(Y5Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_0, MaxAbsValue(Im(Y5Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_0);
   }

   if (MaxAbsValue(Re(Y5Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_1, MaxAbsValue(Re(Y5Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_1);
   }

   if (MaxAbsValue(Im(Y5Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_1, MaxAbsValue(Im(Y5Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_1);
   }

   if (MaxAbsValue(Re(Y5Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_2, MaxAbsValue(Re(Y5Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr0_2);
   }

   if (MaxAbsValue(Im(Y5Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_2, MaxAbsValue(Im(Y5Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr0_2);
   }

   if (MaxAbsValue(Re(Y5Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_0, MaxAbsValue(Re(Y5Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_0);
   }

   if (MaxAbsValue(Im(Y5Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_0, MaxAbsValue(Im(Y5Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_0);
   }

   if (MaxAbsValue(Re(Y5Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_1, MaxAbsValue(Re(Y5Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_1);
   }

   if (MaxAbsValue(Im(Y5Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_1, MaxAbsValue(Im(Y5Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_1);
   }

   if (MaxAbsValue(Re(Y5Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_2, MaxAbsValue(Re(Y5Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr1_2);
   }

   if (MaxAbsValue(Im(Y5Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_2, MaxAbsValue(Im(Y5Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr1_2);
   }

   if (MaxAbsValue(Re(Y5Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_0, MaxAbsValue(Re(Y5Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_0);
   }

   if (MaxAbsValue(Im(Y5Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_0, MaxAbsValue(Im(Y5Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_0);
   }

   if (MaxAbsValue(Re(Y5Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_1, MaxAbsValue(Re(Y5Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_1);
   }

   if (MaxAbsValue(Im(Y5Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_1, MaxAbsValue(Im(Y5Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_1);
   }

   if (MaxAbsValue(Re(Y5Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_2, MaxAbsValue(Re(Y5Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY5Pr2_2);
   }

   if (MaxAbsValue(Im(Y5Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_2, MaxAbsValue(Im(Y5Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY5Pr2_2);
   }

   if (MaxAbsValue(Re(Y1(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10_0, MaxAbsValue(Re(Y1(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10_0);
   }

   if (MaxAbsValue(Im(Y1(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10_0, MaxAbsValue(Im(Y1(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10_0);
   }

   if (MaxAbsValue(Re(Y1(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10_1, MaxAbsValue(Re(Y1(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10_1);
   }

   if (MaxAbsValue(Im(Y1(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10_1, MaxAbsValue(Im(Y1(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10_1);
   }

   if (MaxAbsValue(Re(Y1(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY10_2, MaxAbsValue(Re(Y1(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY10_2);
   }

   if (MaxAbsValue(Im(Y1(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY10_2, MaxAbsValue(Im(Y1(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY10_2);
   }

   if (MaxAbsValue(Re(Y1(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY11_0, MaxAbsValue(Re(Y1(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY11_0);
   }

   if (MaxAbsValue(Im(Y1(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY11_0, MaxAbsValue(Im(Y1(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY11_0);
   }

   if (MaxAbsValue(Re(Y1(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY11_1, MaxAbsValue(Re(Y1(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY11_1);
   }

   if (MaxAbsValue(Im(Y1(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY11_1, MaxAbsValue(Im(Y1(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY11_1);
   }

   if (MaxAbsValue(Re(Y1(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY11_2, MaxAbsValue(Re(Y1(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY11_2);
   }

   if (MaxAbsValue(Im(Y1(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY11_2, MaxAbsValue(Im(Y1(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY11_2);
   }

   if (MaxAbsValue(Re(Y1(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY12_0, MaxAbsValue(Re(Y1(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY12_0);
   }

   if (MaxAbsValue(Im(Y1(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY12_0, MaxAbsValue(Im(Y1(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY12_0);
   }

   if (MaxAbsValue(Re(Y1(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY12_1, MaxAbsValue(Re(Y1(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY12_1);
   }

   if (MaxAbsValue(Im(Y1(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY12_1, MaxAbsValue(Im(Y1(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY12_1);
   }

   if (MaxAbsValue(Re(Y1(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY12_2, MaxAbsValue(Re(Y1(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY12_2);
   }

   if (MaxAbsValue(Im(Y1(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY12_2, MaxAbsValue(Im(Y1(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY12_2);
   }

   if (MaxAbsValue(Re(Y1Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_0, MaxAbsValue(Re(Y1Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_0);
   }

   if (MaxAbsValue(Im(Y1Pr(0,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_0, MaxAbsValue(Im(Y1Pr(0,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_0);
   }

   if (MaxAbsValue(Re(Y1Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_1, MaxAbsValue(Re(Y1Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_1);
   }

   if (MaxAbsValue(Im(Y1Pr(0,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_1, MaxAbsValue(Im(Y1Pr(0,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_1);
   }

   if (MaxAbsValue(Re(Y1Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_2, MaxAbsValue(Re(Y1Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr0_2);
   }

   if (MaxAbsValue(Im(Y1Pr(0,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_2, MaxAbsValue(Im(Y1Pr(0,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr0_2);
   }

   if (MaxAbsValue(Re(Y1Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_0, MaxAbsValue(Re(Y1Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_0);
   }

   if (MaxAbsValue(Im(Y1Pr(1,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_0, MaxAbsValue(Im(Y1Pr(1,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_0);
   }

   if (MaxAbsValue(Re(Y1Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_1, MaxAbsValue(Re(Y1Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_1);
   }

   if (MaxAbsValue(Im(Y1Pr(1,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_1, MaxAbsValue(Im(Y1Pr(1,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_1);
   }

   if (MaxAbsValue(Re(Y1Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_2, MaxAbsValue(Re(Y1Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr1_2);
   }

   if (MaxAbsValue(Im(Y1Pr(1,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_2, MaxAbsValue(Im(Y1Pr(1,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr1_2);
   }

   if (MaxAbsValue(Re(Y1Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_0, MaxAbsValue(Re(Y1Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_0);
   }

   if (MaxAbsValue(Im(Y1Pr(2,0))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_0, MaxAbsValue(Im(Y1Pr(2,0))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_0);
   }

   if (MaxAbsValue(Re(Y1Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_1, MaxAbsValue(Re(Y1Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_1);
   }

   if (MaxAbsValue(Im(Y1Pr(2,1))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_1, MaxAbsValue(Im(Y1Pr(2,1))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_1);
   }

   if (MaxAbsValue(Re(Y1Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_2, MaxAbsValue(Re(Y1Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ReY1Pr2_2);
   }

   if (MaxAbsValue(Im(Y1Pr(2,2))) > perturbativity_bound) {
      problem = true;
      high_scale_model->get_problems().flag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_2, MaxAbsValue(Im(Y1Pr(2,2))),
         high_scale_model->get_scale(), perturbativity_bound);
   } else {
      high_scale_model->get_problems().unflag_non_perturbative_parameter(
         NMFSU5_info::ImY1Pr2_2);
   }

   return problem;
}

double NMFSU5EFT_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double NMFSU5EFT_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const NMFSU5EFT_input_parameters& NMFSU5EFT_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return input;
}

cSMHdCKMRHN<Two_scale>* NMFSU5EFT_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<cSMHdCKMRHN<Two_scale>*>(model_);
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   input = NMFSU5EFT_input_parameters();
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_singlet_yukawas();
   calculate_fiveplet_yukawas();
   calculate_tenplet_yukawas();
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::calculate_singlet_yukawas()
{
   const auto v = HIGHSCALEMODELPARAMETER(v);
   const auto vPr = HIGHSCALEMODELPARAMETER(vPr);
   const auto Y1PrIN = INPUTPARAMETER(Y1PrIN);
   const auto mass_matrix_Fe = model->mass_matrix_Fe();

   Eigen::Matrix<std::complex<double>,3,3> Y1;
   Y1 =( mass_matrix_Fe - Y1PrIN * vPr) / v;

   HIGHSCALEMODEL->set_Y1(Y1);
   HIGHSCALEMODEL->set_Y1Pr(Y1PrIN);
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::calculate_fiveplet_yukawas()
{
   const auto Yu = MODELPARAMETER(Yu);
   MODEL->set_Yv((Yu.transpose()).template cast<std::complex<double> >());

   const auto v = HIGHSCALEMODELPARAMETER(v);
   const auto vPr = HIGHSCALEMODELPARAMETER(vPr);
   const auto Y5bPrIN = INPUTPARAMETER(Y5bPrIN);
   const auto mass_matrix_Fu = model->mass_matrix_Fu();

   Eigen::Matrix<std::complex<double>,3,3> Y5b;
   Y5b = -(Sqrt(2.) * mass_matrix_Fu + Y5bPrIN * vPr) / v;

   HIGHSCALEMODEL->set_Y5b(Y5b);
   HIGHSCALEMODEL->set_Y5bPr(Y5bPrIN);
}

// @todo check signs and prefactors
void NMFSU5EFT_high_scale_constraint<Two_scale>::calculate_tenplet_yukawas()
{
   Eigen::Matrix<std::complex<double>,3,3> Ud;
   Eigen::Matrix<std::complex<double>,3,3> Vd;
   Eigen::Array<double,3,1> Yd_diag;
   fs_svd(MODELPARAMETER(Yd), Yd_diag, Ud, Vd);

   const Eigen::Matrix<std::complex<double>,3,3> Yd_sym(
      Vd.transpose() * Yd_diag.matrix().asDiagonal() * Vd);

   MODEL->set_Yd(Yd_sym);

   const auto mass_matrix_Fd = model->mass_matrix_Fd();
   const auto Mv = MODELPARAMETER(Mv);
   const auto g5 = HIGHSCALEMODELPARAMETER(g5);
   const auto VG = HIGHSCALEMODELPARAMETER(VG);
   const auto MDelta = HIGHSCALEMODELPARAMETER(MDelta);
   const auto MX = HIGHSCALEMODELPARAMETER(MX);
   const auto UDelta = HIGHSCALEMODELPARAMETER(UDelta);
   const auto I3s1 = Sqr(MDelta(0)) / Sqr(MX);
   const auto I3s2 = Sqr(MDelta(1)) / Sqr(MX);
   const auto I3s3 = Sqr(MDelta(2)) / Sqr(MX);

   Eigen::Matrix<std::complex<double>,12,12> lhs(
      Eigen::Matrix<std::complex<double>,12,12>::Zero());
   Eigen::Matrix<std::complex<double>,12,1> rhs(
      Eigen::Matrix<std::complex<double>,12,1>::Zero());

   lhs(0,0) = -0.5 * v;
   lhs(0,6) = -0.5 * vPr;
   lhs(1,1) = -0.5 * v;
   lhs(1,7) = -0.5 * vPr;
   lhs(2,2) = -0.5 * v;
   lhs(2,8) = -0.5 * vPr;
   lhs(3,3) = -0.5 * v;
   lhs(3,9) = -0.5 * vPr;
   lhs(4,4) = -0.5 * v;
   lhs(4,10) = -0.5 * vPr;
   lhs(5,5) = -0.5 * v;
   lhs(5,11) = -0.5 * vPr;
   lhs(6,0) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(6,6) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);
   lhs(7,1) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(7,7) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);
   lhs(8,2) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(8,8) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);
   lhs(9,3) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(9,9) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);
   lhs(10,4) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(10,10) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);
   lhs(11,5) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,1)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,1)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,1)) * I3s3);
   lhs(11,11) = -3. * twoLoop * Pow4(g5) * VG * (
      UDelta(0,0) * Conj(UDelta(0,2)) * I3s1
      + UDelta(1,0) * Conj(UDelta(1,2)) * I3s2
      + UDelta(2,0) * Conj(UDelta(2,2)) * I3s3);

   rhs(0) = mass_matrix_Fd(0,0);
   rhs(1) = mass_matrix_Fd(0,1);
   rhs(2) = mass_matrix_Fd(0,2);
   rhs(3) = mass_matrix_Fd(1,1);
   rhs(4) = mass_matrix_Fd(1,2);
   rhs(5) = mass_matrix_Fd(2,2);
   rhs(6) = Mv(0,0);
   rhs(7) = Mv(0,1);
   rhs(8) = Mv(0,2);
   rhs(9) = Mv(1,1);
   rhs(10) = Mv(1,2);
   rhs(11) = Mv(2,2);

   const Eigen::Matrix<std::complex<double>,12,1> sol
      = lhs.colPivHouseholderQr().solve(rhs);

   Eigen::Matrix<std::complex<double>,3,3> Y10;
   Eigen::Matrix<std::complex<double>,3,3> Y10Pr;

   Y10(0,0) = sol(0);
   Y10(0,1) = sol(1);
   Y10(0,2) = sol(2);
   Y10(1,1) = sol(3);
   Y10(1,2) = sol(4);
   Y10(2,2) = sol(5);
   Y10Pr(0,0) = sol(6);
   Y10Pr(0,1) = sol(7);
   Y10Pr(0,2) = sol(8);
   Y10Pr(1,1) = sol(9);
   Y10Pr(1,2) = sol(10);
   Y10Pr(2,2) = sol(11);

   Symmetrize(Y10);
   Symmetrize(Y10Pr);

   HIGHSCALEMODEL->set_Y10(Y10);
   HIGHSCALEMODEL->set_Y10Pr(Y10Pr);
}

void NMFSU5EFT_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


}

void NMFSU5EFT_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("NMFSU5EFT_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");

   if (!high_scale_model)
      throw SetupError("NMFSU5EFT_high_scale_constraint<Two_scale>: "
                       "high-scale model pointer is zero!");
}

} // namespace flexiblesusy

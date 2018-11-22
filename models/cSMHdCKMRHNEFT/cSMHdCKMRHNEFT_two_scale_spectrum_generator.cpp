#include "cSMHdCKMRHNEFT_two_scale_spectrum_generator.hpp"

#include "cSMHdCKM_input_parameters.hpp"
#include "cSMHdCKM_two_scale_convergence_tester.hpp"
#include "cSMHdCKM_two_scale_model.hpp"

#include "cSMHdCKMRHNEFT_input_parameters.hpp"
#include "cSMHdCKMRHNEFT_two_scale_low_scale_constraint.hpp"
#include "cSMHdCKMRHNEFT_two_scale_matching.hpp"
#include "cSMHdCKMRHNEFT_two_scale_susy_scale_constraint.hpp"

#include "cSMHdCKMRHN_two_scale_convergence_tester.hpp"
#include "cSMHdCKMRHN_two_scale_ewsb_solver.hpp"
#include "cSMHdCKMRHNEFT_two_scale_high_scale_constraint.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"

#include "cSMHdCKMRHNEFT_two_scale_initial_guesser.hpp"

#include "composite_convergence_tester.hpp"
#include "error.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

namespace flexiblesusy {

double cSMHdCKMRHNEFT_spectrum_generator<Two_scale>::get_pole_mass_scale(double susy_scale) const
{
   return settings.get(Spectrum_generator_settings::pole_mass_scale) != 0. ?
      settings.get(Spectrum_generator_settings::pole_mass_scale) :
      susy_scale;
}

double cSMHdCKMRHNEFT_spectrum_generator<Two_scale>::get_eft_pole_mass_scale(double susy_scale, double Mt) const
{
   double Q_higgs = settings.get(
      Spectrum_generator_settings::eft_pole_mass_scale);

   if (Q_higgs == 0.)
      Q_higgs = std::min(susy_scale, Mt);

   return Q_higgs;
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale,
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
void cSMHdCKMRHNEFT_spectrum_generator<Two_scale>::run_except(
   const softsusy::QedQcd& qedqcd,
   const cSMHdCKMRHNEFT_input_parameters& input)
{
   VERBOSE_MSG("Solving BVP using two-scale solver");

   problems.set_bvp_solver_problems({ BVP_solver_problems("TwoScaleSolver") });

   auto& model = this->model;
   model.clear();
   model.do_calculate_sm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_sm_masses));
   model.do_calculate_bsm_pole_masses(
      settings.get(Spectrum_generator_settings::calculate_bsm_masses));
   model.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   model.set_loops(settings.get(Spectrum_generator_settings::beta_loop_order));
   model.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   model.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));

   eft.clear();
   eft.do_calculate_sm_pole_masses(true);
   eft.do_force_output(
      settings.get(Spectrum_generator_settings::force_output));
   eft.set_loops(settings.get(Spectrum_generator_settings::beta_loop_order));
   eft.set_thresholds(
      settings.get(
         Spectrum_generator_settings::threshold_corrections_loop_order));
   eft.set_zero_threshold(
      settings.get(Spectrum_generator_settings::beta_zero_threshold));
   eft.set_pole_mass_loop_order(this->model.get_pole_mass_loop_order());
   eft.set_ewsb_loop_order(this->model.get_ewsb_loop_order());
   eft.set_ewsb_iteration_precision(this->model.get_ewsb_iteration_precision());
   eft.set_loop_corrections(this->model.get_loop_corrections());
   eft.set_threshold_corrections(this->model.get_threshold_corrections());

   cSMHdCKMRHN_ewsb_solver<Two_scale> ewsb_solver;
   model.set_ewsb_solver(
      std::make_shared<cSMHdCKMRHN_ewsb_solver<Two_scale> >(ewsb_solver));

   cSMHdCKMRHNEFT_high_scale_constraint<Two_scale> high_scale_constraint(&model, input);
   cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale> susy_scale_constraint(&model, qedqcd, input);
   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&eft, qedqcd, input);

   // note: to avoid large logarithms the downwards matching loop order
   // is used for both matching conditions
   const int matching_loop_order_up = settings.get(
      Spectrum_generator_settings::eft_matching_loop_order_down);
   const int matching_loop_order_down = settings.get(
      Spectrum_generator_settings::eft_matching_loop_order_down);
   const auto scale_getter = [this,&susy_scale_constraint] () {
      return susy_scale_constraint.get_scale(); };

   cSMHdCKMRHNEFT_matching_up<Two_scale> matching_up(
      &eft, &model, scale_getter, matching_loop_order_up);
   cSMHdCKMRHNEFT_matching_down<Two_scale> matching_down(
      &eft, &model, scale_getter, matching_loop_order_down);

   matching_up.set_scale(
      settings.get(Spectrum_generator_settings::eft_matching_scale));
   matching_down.set_scale(
      settings.get(Spectrum_generator_settings::eft_matching_scale));

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();

   // convergence tester for cSMHdCKMRHN
   cSMHdCKMRHN_convergence_tester<Two_scale> model_ct(
      &model, settings.get(Spectrum_generator_settings::precision),
      [this,&susy_scale_constraint](){
         return get_pole_mass_scale(susy_scale_constraint.get_scale()); });

   // convergence tester for cSMHdCKM
   const double Mt = qedqcd.displayPoleMt();
   cSMHdCKM_convergence_tester<Two_scale> eft_ct(
      &eft, settings.get(Spectrum_generator_settings::precision),
      [this,&susy_scale_constraint,Mt](){
         return get_eft_pole_mass_scale(susy_scale_constraint.get_scale(), Mt);
      });

   if (settings.get(Spectrum_generator_settings::max_iterations) > 0) {
      model_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
      eft_ct.set_max_iterations(
         settings.get(Spectrum_generator_settings::max_iterations));
   }

   Composite_convergence_tester cct;
   cct.add_convergence_tester(&model_ct);
   cct.add_convergence_tester(&eft_ct);

   cSMHdCKMRHNEFT_initial_guesser<Two_scale> initial_guesser(
      &model, &eft, qedqcd,
      low_scale_constraint, susy_scale_constraint, high_scale_constraint);

   Two_scale_increasing_precision precision(
      10.0, settings.get(Spectrum_generator_settings::precision));

   RGFlow<Two_scale> solver;
   solver.set_convergence_tester(&cct);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add(&low_scale_constraint, &eft);
   solver.add(&matching_up, &eft, &model);
   solver.add(&susy_scale_constraint, &model);
   solver.add(&high_scale_constraint, &model);
   solver.add(&susy_scale_constraint, &model);
   solver.add(&matching_down, &model, &eft);

   high_scale = susy_scale = low_scale = 0.;
   reached_precision = std::numeric_limits<double>::infinity();

   solver.solve();

   // impose low-scale constraint one last time
   eft.run_to(low_scale_constraint.get_scale());
   low_scale_constraint.apply();

   high_scale = high_scale_constraint.get_scale();
   susy_scale = susy_scale_constraint.get_scale();
   low_scale  = low_scale_constraint.get_scale();
   reached_precision = std::max(model_ct.get_current_accuracy(),
                                eft_ct.get_current_accuracy());

   calculate_spectrum(Mt, low_scale_constraint.get_sm_parameters().displayPoleMW());

   // run to output scale (if scale > 0)
   if (!is_zero(parameter_output_scale)) {
      model.run_to(parameter_output_scale);
      eft.run_to(parameter_output_scale);
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
void cSMHdCKMRHNEFT_spectrum_generator<Two_scale>::write_running_couplings(
   const std::string& filename) const
{
   cSMHdCKMRHNEFT_spectrum_generator_interface<Two_scale>::write_running_couplings(
      filename, get_low_scale(), get_high_scale());
}

void cSMHdCKMRHNEFT_spectrum_generator<Two_scale>::calculate_spectrum(double Mt, double MW)
{
   model.run_to(get_pole_mass_scale(susy_scale));
   model.solve_ewsb();
   model.calculate_spectrum();

   eft.run_to(get_eft_pole_mass_scale(susy_scale, Mt));

   // computation of pole mass spectrum in the SM
   const int eft_pole_loops = eft.get_pole_mass_loop_order();
   const int eft_ewsb_loops = eft.get_ewsb_loop_order();

   eft.calculate_DRbar_masses();
   eft.solve_ewsb();
   eft.calculate_spectrum();

   eft.set_pole_mass_loop_order(eft_pole_loops);
   eft.set_ewsb_loop_order(eft_ewsb_loops);

   model.get_physical().Mhh = eft.get_physical().Mhh;
   model.get_physical().MVZ = eft.get_physical().MVZ;
   model.get_physical().MVWp = MW;
   model.get_physical().MFu = eft.get_physical().MFu;
   model.get_physical().MFd = eft.get_physical().MFd;
   model.get_physical().MFe = eft.get_physical().MFe;
   model.get_physical().MFv(0) = eft.get_physical().MFv(0);
   model.get_physical().MFv(1) = eft.get_physical().MFv(1);
   model.get_physical().MFv(2) = eft.get_physical().MFv(2);

   if (eft.get_problems().is_running_tachyon(cSMHdCKM_info::hh))
      model.get_problems().flag_running_tachyon(cSMHdCKMRHN_info::hh);
   if (eft.get_problems().is_pole_tachyon(cSMHdCKM_info::hh))
      model.get_problems().flag_pole_tachyon(cSMHdCKMRHN_info::hh);
   if (eft.get_problems().is_running_tachyon(cSMHdCKM_info::VZ))
      model.get_problems().flag_running_tachyon(cSMHdCKMRHN_info::VZ);
   if (eft.get_problems().is_pole_tachyon(cSMHdCKM_info::VZ))
      model.get_problems().flag_pole_tachyon(cSMHdCKMRHN_info::VZ);
   if (eft.get_problems().is_running_tachyon(cSMHdCKM_info::VWp))
      model.get_problems().flag_running_tachyon(cSMHdCKMRHN_info::VWp);
   if (eft.get_problems().is_pole_tachyon(cSMHdCKM_info::VWp))
      model.get_problems().flag_pole_tachyon(cSMHdCKMRHN_info::VWp);
}

} // namespace flexiblesusy

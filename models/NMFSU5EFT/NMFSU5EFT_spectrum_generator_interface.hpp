#ifndef NMFSU5EFT_SPECTRUM_GENERATOR_INTERFACE_H
#define NMFSU5EFT_SPECTRUM_GENERATOR_INTERFACE_H

#include "cSMHdCKM_mass_eigenstates.hpp"
#include "cSMHdCKM_model.hpp"
#include "cSMHdCKM_model_slha.hpp"

#include "cSMHdCKMRHN_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_model.hpp"
#include "cSMHdCKMRHN_model_slha.hpp"
#include "cSMHdCKMRHN_utilities.hpp"

#include "NMFSU5_mass_eigenstates.hpp"
#include "NMFSU5_model.hpp"
#include "NMFSU5_model_slha.hpp"

#include "coupling_monitor.hpp"
#include "error.hpp"
#include "lowe.h"
#include "spectrum_generator_problems.hpp"
#include "spectrum_generator_settings.hpp"
#include "loop_corrections.hpp"

#include <tuple>

namespace flexiblesusy {

struct NMFSU5EFT_input_parameters;

template <class T>
class NMFSU5EFT_spectrum_generator_interface {
public:
   virtual ~NMFSU5EFT_spectrum_generator_interface() = default;

   std::tuple<NMFSU5<T>, cSMHdCKMRHN<T>,cSMHdCKM<T> > get_models() const {
      return std::make_tuple(high_scale_model, model, eft);
   }
   std::tuple<
      NMFSU5_slha<NMFSU5<T> >,
      cSMHdCKMRHN_slha<cSMHdCKMRHN<T> >,
      cSMHdCKM_slha<cSMHdCKM<T> > >
   get_models_slha() const {
      return std::make_tuple(
         NMFSU5_slha<NMFSU5<T> >(high_scale_model, settings.get(
                                    Spectrum_generator_settings::force_positive_masses) == 0.),
         cSMHdCKMRHN_slha<cSMHdCKMRHN<T> >(model, settings.get(
                                        Spectrum_generator_settings::force_positive_masses) == 0.),
         cSMHdCKM_slha<cSMHdCKM<T> >(eft, settings.get(
                                        Spectrum_generator_settings::force_positive_masses) == 0.));
   }

   const NMFSU5<T>& get_high_scale_model() const { return high_scale_model; }
   NMFSU5<T>& get_high_scale_model() const { return high_scale_model; }
   NMFSU5_slha<NMFSU5<T> > get_high_scale_model_slha() const {
      return NMFSU5_slha<NMFSU5<T> >(
         high_scale_model, settings.get(
            Spectrum_generator_settings::force_positive_masses) == 0.);
   }

   const cSMHdCKMRHN<T>& get_model() const { return model; }
   cSMHdCKMRHN<T>& get_model() { return model; }
   cSMHdCKMRHN_slha<cSMHdCKMRHN<T> > get_model_slha() const {
      return cSMHdCKMRHN_slha<cSMHdCKMRHN<T> >(
         model, settings.get(
            Spectrum_generator_settings::force_positive_masses) == 0.);
   }

   const cSMHdCKM<T>& get_eft() const { return eft; }
   cSMHdCKM<T>& get_eft() { return eft; }
   cSMHdCKM_slha<cSMHdCKM<T> > get_eft_slha() const {
      return cSMHdCKM_slha<cSMHdCKM<T> >(
         eft, settings.get(Spectrum_generator_settings::force_positive_masses) == 0.);
   }

   Spectrum_generator_problems get_problems() const { return problems; }
   int get_exit_code() const { return problems.have_problem(); }
   double get_reached_precision() const { return reached_precision; }
   const Spectrum_generator_settings& get_settings() const { return settings; }
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_settings(const Spectrum_generator_settings&);

   void run(const softsusy::QedQcd&, const NMFSU5EFT_input_parameters&);
   void write_running_couplings(const std::string& filename, double, double) const;
   void write_spectrum(const std::string& filename = "NMFSU5EFT_spectrum.dat") const;

protected:
   NMFSU5<T> high_scale_model;
   cSMHdCKMRHN<T> model;
   cSMHdCKM<T> eft;
   Spectrum_generator_problems problems;
   Spectrum_generator_settings settings;
   double parameter_output_scale{0.}; ///< output scale for running parameters
   double reached_precision{std::numeric_limits<double>::infinity()}; ///< the precision that was reached

   void translate_exception_to_problem(NMFSU5<T>&);
   virtual void run_except(
      const softsusy::QedQcd&, const NMFSU5EFT_input_parameters&) = 0;
};

/**
 * Setup spectrum generator from a Spectrum_generator_settings object.
 *
 * @param settings_ spectrum generator settings
 */
template <class T>
void NMFSU5EFT_spectrum_generator_interface<T>::set_settings(
   const Spectrum_generator_settings& settings_)
{
   settings = settings_;

   high_scale_model.set_pole_mass_loop_order(
      settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   high_scale_model.set_ewsb_loop_order(
      settings.get(Spectrum_generator_settings::ewsb_loop_order));
   high_scale_model.set_loop_corrections(settings.get_loop_corrections());
   high_scale_model.set_threshold_corrections(
      settings.get_threshold_corrections());

   model.set_pole_mass_loop_order(settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   model.set_ewsb_loop_order(settings.get(Spectrum_generator_settings::ewsb_loop_order));
   model.set_loop_corrections(settings.get_loop_corrections());
   model.set_threshold_corrections(settings.get_threshold_corrections());

   eft.set_pole_mass_loop_order(settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   eft.set_ewsb_loop_order(settings.get(Spectrum_generator_settings::ewsb_loop_order));
   eft.set_loop_corrections(settings.get_loop_corrections());
   eft.set_threshold_corrections(settings.get_threshold_corrections());
}

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function calls run_except() from the derived class and
 * translates an emitted exception into an problem code.
 *
 * @param qedqcd_ Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void NMFSU5EFT_spectrum_generator_interface<T>::run(
   const softsusy::QedQcd& qedqcd_, const NMFSU5EFT_input_parameters& input)
{
   softsusy::QedQcd qedqcd = qedqcd_;

   try {
      qedqcd.to(qedqcd.displayPoleMZ());
      this->run_except(qedqcd, input);
   } catch (...) {
      this->translate_exception_to_problem(model);
   }

   problems.set_model_problems(
      { high_scale_model.get_problems(), model.get_problems(), eft.get_problems() });
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 * @param start lowest scale
 * @param stop highest scale
 */
template <class T>
void NMFSU5EFT_spectrum_generator_interface<T>::write_running_couplings(
   const std::string& filename,
   double start, double stop) const
{
   NMFSU5_mass_eigenstates tmp_model(high_scale_model);
   try {
      tmp_model.run_to(start);
   } catch (const Error& error) {
      ERROR("write_running_couplings: running to scale "
            << start << " failed: " << error.what());
      return;
   }

   NMFSU5_parameter_getter parameter_getter;
   Coupling_monitor<NMFSU5_mass_eigenstates, NMFSU5_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(start, stop, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void NMFSU5EFT_spectrum_generator_interface<T>::write_spectrum(
   const std::string& filename) const
{
   NMFSU5_spectrum_plotter plotter(high_scale_model);
   plotter.write_to_file(filename);
}

/**
 * Flags problems in the given model class from the current pending
 * exception.
 *
 * This function assumes that there is an active exception.
 *
 * @param model model class
 */
template <class T>
void NMFSU5EFT_spectrum_generator_interface<T>::translate_exception_to_problem(
   NMFSU5<T>& model)
{
   try {
      throw;
   } catch (const NoConvergenceError&) {
      problems.flag_no_convergence();
   } catch (const NonPerturbativeRunningError& error) {
      high_scale_model.get_problems().flag_no_perturbative();
      model.get_problems().flag_no_perturbative();
      eft.get_problems().flag_no_perturbative();
      high_scale_model.get_problems().flag_non_perturbative_parameter(
         error.get_parameter_index(), error.get_parameter_value(),
         error.get_scale());
   } catch (const NonPerturbativeRunningQedQcdError& error) {
      high_scale_model.get_problems().flag_no_perturbative();
      high_scale_model.get_problems().flag_thrown(error.what());
      model.get_problems().flag_no_perturbative();
      model.get_problems().flag_thrown(error.what());
      eft.get_problems().flag_no_perturbative();
      eft.get_problems().flag_thrown(error.what());
   } catch (const NoSinThetaWConvergenceError&) {
      high_scale_model.get_problems().flag_no_sinThetaW_convergence();
      model.get_problems().flag_no_sinThetaW_convergence();
      eft.get_problems().flag_no_sinThetaW_convergence();
   } catch (const SetupError& error) {
      high_scale_model.get_problems().flag_thrown(error.what());
      model.get_problems().flag_thrown(error.what());
      eft.get_problems().flag_thrown(error.what());
   } catch (const Error& error) {
      high_scale_model.get_problems().flag_thrown(error.what());
      model.get_problems().flag_thrown(error.what());
      eft.get_problems().flag_thrown(error.what());
   } catch (const std::exception& error) {
      high_scale_model.get_problems().flag_thrown(error.what());
      model.get_problems().flag_thrown(error.what());
      eft.get_problems().flag_thrown(error.what());
   }
}

} // namespace flexiblesusy

#endif

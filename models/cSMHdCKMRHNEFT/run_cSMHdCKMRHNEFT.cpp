
#include "config.h"

#include "cSMHdCKM_input_parameters.hpp"

#include "cSMHdCKMRHN_input_parameters.hpp"
#include "cSMHdCKMRHN_observables.hpp"
#include "cSMHdCKMRHN_utilities.hpp"

#include "cSMHdCKMRHNEFT_slha_io.hpp"
#include "cSMHdCKMRHNEFT_spectrum_generator.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "cSMHdCKMRHNEFT_two_scale_spectrum_generator.hpp"
#endif

#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

/**
 * @brief Runs the spectrum generator of type \a solver_type
 * @tparam solver_type solver type
 * @param slha_io SLHA input
 * @param spectrum_generator_settings
 * @param slha_output_file output file for SLHA output
 * @param database_output_file output file for SQLite database
 * @param spectrum_file output file for the mass spectrum
 * @param rgflow_file output file for the RG flow
 * @return value of spectrum_generator::get_exit_code()
 */
template <class solver_type>
int run_solver(flexiblesusy::cSMHdCKMRHNEFT_slha_io& slha_io,
               const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings,
               const std::string& slha_output_file,
               const std::string& database_output_file,
               const std::string& spectrum_file,
               const std::string& rgflow_file)
{
   using namespace flexiblesusy;

   Physical_input physical_input; // extra non-SLHA physical input
   softsusy::QedQcd qedqcd;
   cSMHdCKM_input_parameters eft_input;
   cSMHdCKMRHN_input_parameters model_input;

   try {
      slha_io.fill(qedqcd);
      slha_io.fill(eft_input, model_input);
      slha_io.fill(physical_input);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   cSMHdCKMRHNEFT_spectrum_generator<solver_type> spectrum_generator;
   spectrum_generator.set_settings(spectrum_generator_settings);
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());

   spectrum_generator.run(qedqcd, eft_input, model_input);

   auto models = spectrum_generator.get_models_slha();
   const auto& problems = spectrum_generator.get_problems();

   cSMHdCKMRHNEFT_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale  = spectrum_generator.get_low_scale();
   scales.pole_mass_scale = spectrum_generator.get_pole_mass_scale();

   cSMHdCKMRHN_observables observables;
   if (spectrum_generator_settings.get(Spectrum_generator_settings::calculate_observables))
      observables = calculate_observables(std::get<0>(models), qedqcd, physical_input, scales.pole_mass_scale);

   const bool show_result = !problems.have_problem() ||
      spectrum_generator_settings.get(Spectrum_generator_settings::force_output);
   // SLHA output
   if (!slha_output_file.empty()) {
      slha_io.set_spinfo(problems);
      slha_io.set_input(eft_input, model_input);
      if (show_result) {
         slha_io.set_print_imaginary_parts_of_majorana_mixings(
            spectrum_generator_settings.get(
               Spectrum_generator_settings::force_positive_masses));
         slha_io.set_spectrum(models);
         slha_io.set_extra(std::get<0>(models), scales, observables);
      }

      slha_io.write_to(slha_output_file);
   }

   if (!database_output_file.empty() && show_result) {
      cSMHdCKMRHN_database::to_database(
         database_output_file, std::get<0>(models), &qedqcd,
         &physical_input, &observables);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   return spectrum_generator.get_exit_code();
}

/**
 * @brief Runs the spectrum generator
 *
 * Reads the solver type from \a spectrum_generator_settings and calls
 * run_solver() with the corresponding solver type.
 *
 * @param slha_io SLHA input
 * @param spectrum_generator_settings
 * @param slha_output_file output file for SLHA output
 * @param database_output_file output file for SQLite database
 * @param spectrum_file output file for the mass spectrum
 * @param rgflow_file output file for the RG flow
 * @return return value of run_solver<>()
 */
int run(
   flexiblesusy::cSMHdCKMRHNEFT_slha_io& slha_io,
   const flexiblesusy::Spectrum_generator_settings& spectrum_generator_settings,
   const std::string& slha_output_file,
   const std::string& database_output_file,
   const std::string& spectrum_file,
   const std::string& rgflow_file)
{
   using namespace flexiblesusy;

   int exit_code = 0;
   const int solver_type
      = static_cast<int>(spectrum_generator_settings.get(
                            Spectrum_generator_settings::solver));

   switch (solver_type) {
   case 0:
#ifdef ENABLE_TWO_SCALE_SOLVER
   case 1:
      exit_code = run_solver<Two_scale>(
         slha_io, spectrum_generator_settings, slha_output_file,
         database_output_file, spectrum_file, rgflow_file);
      if (!exit_code || solver_type != 0) break;
#endif

   default:
      if (solver_type != 0) {
         ERROR("unknown solver type: " << solver_type);
         exit_code = -1;
      }
      break;
   }

   return exit_code;
}

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      cSMHdCKMRHN_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string database_output_file(options.get_database_output_file());
   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_source(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   cSMHdCKMRHNEFT_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;

   if (slha_input_source.empty()) {
      ERROR("No SLHA input source given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_source(slha_input_source);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   const int exit_code
      = run(slha_io, spectrum_generator_settings, slha_output_file,
            database_output_file, spectrum_file, rgflow_file);

   return exit_code;
}


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cSMHdCKMRHNEFT_boundary_conditions

#include "cSMHdCKMRHNEFT_input_parameters.hpp"
#include "cSMHdCKM_model_slha.hpp"
#include "cSMHdCKMRHN_model_slha.hpp"
#include "cSMHdCKMRHNEFT_two_scale_spectrum_generator.hpp"

#include "matrix_tests.hpp"
#include "random_matrix.hpp"

#include <boost/test/unit_test.hpp>

#include <tuple>

using namespace flexiblesusy;

void initialize_sm_input(softsusy::QedQcd& qedqcd)
{
   CKM_parameters ckm;
   ckm.reset_to_observation();

   PMNS_parameters pmns;
   pmns.reset_to_observation();

   qedqcd.setCKM(ckm);
   qedqcd.setPMNS(pmns);

   qedqcd.setNeutrinoPoleMass(1, 1.e-12);
   qedqcd.setNeutrinoPoleMass(2, 5.e-12);
   qedqcd.setNeutrinoPoleMass(3, 1.e-11);
}

void initialize_model_input(cSMHdCKMRHNEFT_input_parameters& input)
{
   input.LambdaIN = 0.2;
   input.Qin = 2.e16;
   input.sign_delta_mAsq = 1;

#ifdef ENABLE_RANDOM
   std::mt19937 generator;
   random_cue_matrix(input.UvInput, generator);
#endif
}

struct Spectrum_generator_result {
   std::tuple<cSMHdCKMRHN_slha<cSMHdCKMRHN<Two_scale> >, cSMHdCKM_slha<cSMHdCKM<Two_scale> > > models;
   Spectrum_generator_problems problems;
   double high_scale;
   double low_scale;
};

Spectrum_generator_result run_spectrum_generator(
   const softsusy::QedQcd& qedqcd, const cSMHdCKMRHNEFT_input_parameters& input)
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 0);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 0);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
   settings.set(Spectrum_generator_settings::precision, 1.0e-4);

   cSMHdCKMRHNEFT_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   return { spectrum_generator.get_models_slha(), spectrum_generator.get_problems(),
         spectrum_generator.get_high_scale(), spectrum_generator.get_low_scale()};
}

BOOST_AUTO_TEST_CASE( test_high_scale_Yd_symmetric )
{
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;
   initialize_sm_input(qedqcd);
   initialize_model_input(input);

   const auto result = run_spectrum_generator(qedqcd, input);

   const bool error = result.problems.have_problem();
   BOOST_REQUIRE(!error);

   if (!error) {
      auto high_scale_model = std::get<0>(result.models);
      const double high_scale = result.high_scale;

      high_scale_model.run_to(high_scale);

      BOOST_CHECK(is_symmetric(high_scale_model.get_Yd(), 1.e-10));
   }
}

BOOST_AUTO_TEST_CASE( test_high_scale_Yu_equal_Yv_transpose )
{
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;
   initialize_sm_input(qedqcd);
   initialize_model_input(input);

   const auto result = run_spectrum_generator(qedqcd, input);

   const bool error = result.problems.have_problem();
   BOOST_REQUIRE(!error);

   if (!error) {
      auto high_scale_model = std::get<0>(result.models);
      const double high_scale = result.high_scale;

      high_scale_model.run_to(high_scale);

      BOOST_CHECK(is_equal(high_scale_model.get_Yv(), high_scale_model.get_Yu().transpose(), 1.e-10));
   }
}

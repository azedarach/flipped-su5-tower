
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cSMHdCKMRHNEFT_low_scale_constraint

#include "matrix_tests.hpp"

#include <boost/test/unit_test.hpp>

#include <cstdlib>

#include "cSMHdCKM_two_scale_model.hpp"

#define private public

#include "cSMHdCKMRHNEFT_two_scale_low_scale_constraint.hpp"

using namespace flexiblesusy;

double random_angle()
{
   return 2. * M_PI * rand() / RAND_MAX;
}

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
   input.UV_theta21 = random_angle();
   input.UV_theta31 = random_angle();
   input.UV_theta32 = random_angle();
   input.UV_phi21 = random_angle();
   input.UV_phi31 = random_angle();
   input.UV_phi32 = random_angle();
   input.UV_chi21 = random_angle();
   input.UV_chi32 = random_angle();
   input.UV_gamma = random_angle();
}

void initialize_model(cSMHdCKM<Two_scale>& model)
{

}

// test that construction of neutrino mixing matrix from
// input angles is unitary
BOOST_AUTO_TEST_CASE( test_neutrino_mixing_matrix_construction )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);
   low_scale_constraint.initialize();

   low_scale_constraint.calculate_neutrino_mixings();

   BOOST_CHECK(is_unitary(low_scale_constraint.neutrinoMix, 1.e-12));
}

// test that the constructed up-type quark Yukawas are
// diagonal and that the resulting CKM mixing matches
// the input CKM matrix up to changes in phase convention
BOOST_AUTO_TEST_CASE( test_quark_mixing )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);
   low_scale_constraint.initialize();

   low_scale_constraint.apply();

   BOOST_CHECK(is_diagonal(model.get_Yu()));

   Eigen::Matrix<std::complex<double>,3,3> Vu(model.get_Vu());
   Eigen::Matrix<std::complex<double>,3,3> Uu(model.get_Uu());
   Eigen::Matrix<std::complex<double>,3,3> Vd(model.get_Vd());
   Eigen::Matrix<std::complex<double>,3,3> Ud(model.get_Ud());
   Eigen::Matrix<std::complex<double>,3,3> ckm(Vu * Vd.adjoint());

   CKM_parameters::to_pdg_convention(ckm, Vu, Vd, Uu, Ud);

   BOOST_CHECK(is_equal(ckm, low_scale_constraint.get_ckm(), 1.e-10));
}

// test that the constructed PMNS matrix agrees with
// the input value up to changes in phase convention
// in the case of normal hierarchy
BOOST_AUTO_TEST_CASE( test_lepton_mixing_nh )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

}

// test that the constructed PMNS matrix agrees with
// the input value up to changes in phase convention
// in the case of inverse hierarchy
BOOST_AUTO_TEST_CASE( test_lepton_mixing_ih )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

}


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cSMHdCKMRHNEFT_low_scale_constraint

#include "config.h"
#include "ew_input.hpp"

#include "matrix_tests.hpp"
#include "random_matrix.hpp"

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#ifdef ENABLE_RANDOM
#include <random>
#endif

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

#ifdef ENABLE_RANDOM
   std::mt19937 generator;
   random_cue_matrix(input.UvInput, generator);
#endif
}

void initialize_model(cSMHdCKM<Two_scale>& model)
{
   model.set_scale(Electroweak_constants::MZ);
   model.set_loops(2);

   model.set_g1(Electroweak_constants::g1);
   model.set_g2(Electroweak_constants::g2);
   model.set_g3(Electroweak_constants::g3);

   model.set_Yu(0, 0, Electroweak_constants::yuSM);
   model.set_Yu(1, 1, Electroweak_constants::ycSM);
   model.set_Yu(2, 2, Electroweak_constants::ytSM);

   model.set_Yd(0, 0, Electroweak_constants::ydSM);
   model.set_Yd(1, 1, Electroweak_constants::ysSM);
   model.set_Yd(2, 2, Electroweak_constants::ybSM);

   model.set_Ye(0, 0, Electroweak_constants::yeSM);
   model.set_Ye(1, 1, Electroweak_constants::ymSM);
   model.set_Ye(2, 2, Electroweak_constants::ylSM);

   model.set_Lambdax(Electroweak_constants::lamSM);

   model.set_mu2(Electroweak_constants::mu2SM);
   model.set_v(Electroweak_constants::vev);
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

   model.calculate_DRbar_masses();

   BOOST_CHECK(is_diagonal(model.get_Yu()));

   Eigen::Matrix<std::complex<double>,3,3> Vu(model.get_Vu());
   Eigen::Matrix<std::complex<double>,3,3> Uu(model.get_Uu());
   Eigen::Matrix<std::complex<double>,3,3> Vd(model.get_Vd());
   Eigen::Matrix<std::complex<double>,3,3> Ud(model.get_Ud());
   Eigen::Matrix<std::complex<double>,3,3> ckm(Vu * Vd.adjoint());

   CKM_parameters::to_pdg_convention(ckm, Vu, Vd, Uu, Ud);

   BOOST_CHECK(is_equal(ckm, low_scale_constraint.get_ckm(), 1.e-10));
}

// test that the input Uv matrix diagonalises the calculated
// neutrino coupling matrix
BOOST_AUTO_TEST_CASE( test_Kappa_diagonalization )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

   model.calculate_DRbar_masses();

   const Eigen::Matrix<std::complex<double>,3,3>& Kappa = model.get_Kappa();
   const Eigen::Matrix<std::complex<double>,3,3>& Kappa_diag
      = input.UvInput.conjugate() * Kappa * input.UvInput.adjoint();

   // scale to obtain non-negligible matrix entries
   const double scaling = 1.e10;
   const Eigen::Matrix<std::complex<double>, 3, 3> scaled_Kappa_diag = scaling * Kappa_diag;

   BOOST_CHECK(is_diagonal(scaled_Kappa_diag, 1.e-15));
}

// tests that in the case of normal hierarchy, the change of neutrino
// basis is trivial
BOOST_AUTO_TEST_CASE( test_nh_basis_change )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

   model.calculate_DRbar_masses();

   const auto mass_matrix_Fv = model.get_mass_matrix_Fv();
   const auto MFv = model.get_MFv();

   const Eigen::Matrix<double,3,3> mass_matrix_diag
      = (low_scale_constraint.neutrinoMix.conjugate() * mass_matrix_Fv
         * low_scale_constraint.neutrinoMix.adjoint()).real();

   BOOST_CHECK(is_equal(mass_matrix_diag(0,0), MFv(0), 1.e-10));
   BOOST_CHECK(is_equal(mass_matrix_diag(1,1), MFv(1), 1.e-10));
   BOOST_CHECK(is_equal(mass_matrix_diag(2,2), MFv(2), 1.e-10));
}

// tests that in the case of inverted hierarchy, the conversion from the
// mass ordered basis to the traditional PMNS basis is performed correctly
BOOST_AUTO_TEST_CASE( test_ih_basis_change )
{
   cSMHdCKM<Two_scale> model;
   softsusy::QedQcd qedqcd;
   cSMHdCKMRHNEFT_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   input.sign_delta_mAsq = -1;

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

   model.calculate_DRbar_masses();

   const auto mass_matrix_Fv = model.get_mass_matrix_Fv();
   const auto MFv = model.get_MFv();

   const Eigen::Matrix<double,3,3> mass_matrix_diag
      = (low_scale_constraint.neutrinoMix.conjugate() * mass_matrix_Fv
         * low_scale_constraint.neutrinoMix.adjoint()).real();

   BOOST_CHECK(is_equal(mass_matrix_diag(0,0), MFv(1), 1.e-10));
   BOOST_CHECK(is_equal(mass_matrix_diag(1,1), MFv(2), 1.e-10));
   BOOST_CHECK(is_equal(mass_matrix_diag(2,2), MFv(0), 1.e-10));
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

   model.calculate_DRbar_masses();

   Eigen::Matrix<std::complex<double>,3,3> Ve(model.get_Ve());
   Eigen::Matrix<std::complex<double>,3,3> Ue(model.get_Ue());
   Eigen::Matrix<std::complex<double>,3,3> UV(model.get_UV());
   Eigen::Matrix<std::complex<double>,3,3> pmns(Ve * UV.adjoint());
   const Eigen::Matrix<double,3,3> pmns_squared_invariants(pmns.cwiseAbs2());

   BOOST_CHECK(is_unitary(UV, 1.e-12));

   PMNS_parameters::to_pdg_convention(pmns, UV, Ve, Ue);

   const Eigen::Matrix<double,3,3> pdg_squared_invariants
      = pmns.cwiseAbs2();

   BOOST_CHECK(is_equal(pmns_squared_invariants, pdg_squared_invariants, 1.e-10));

   const Eigen::Matrix<std::complex<double>,3,3> expected_pmns
      = low_scale_constraint.get_pmns();
   const Eigen::Matrix<double,3,3> expected_squared_invariants
      = expected_pmns.cwiseAbs2();

   BOOST_CHECK(is_equal(pmns, expected_pmns, 1.e-10));
   BOOST_CHECK(is_equal(pmns_squared_invariants, expected_squared_invariants, 1.e-10));
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

   input.sign_delta_mAsq = -1;

   cSMHdCKMRHNEFT_low_scale_constraint<Two_scale> low_scale_constraint(&model, qedqcd, input);

   low_scale_constraint.apply();

   model.calculate_DRbar_masses();

   Eigen::Matrix<std::complex<double>,3,3> Ve(model.get_Ve());
   Eigen::Matrix<std::complex<double>,3,3> Ue(model.get_Ue());
   // note in this case the value of UV from the model cannot be
   // used, as the basis is always assumed to be mass ordered
   Eigen::Matrix<std::complex<double>,3,3> UV(low_scale_constraint.neutrinoMix);
   Eigen::Matrix<std::complex<double>,3,3> pmns(Ve * UV.adjoint());

   const Eigen::Matrix<double,3,3> pmns_squared_invariants(pmns.cwiseAbs2());

   BOOST_CHECK(is_unitary(UV, 1.e-12));

   PMNS_parameters::to_pdg_convention(pmns, UV, Ve, Ue);

   const Eigen::Matrix<double,3,3> pdg_squared_invariants
      = pmns.cwiseAbs2();

   BOOST_CHECK(is_equal(pmns_squared_invariants, pdg_squared_invariants, 1.e-10));

   const Eigen::Matrix<std::complex<double>,3,3> expected_pmns
      = low_scale_constraint.get_pmns();
   const Eigen::Matrix<double,3,3> expected_squared_invariants
      = expected_pmns.cwiseAbs2();

   BOOST_CHECK(is_equal(pmns, expected_pmns, 1.e-10));
   BOOST_CHECK(is_equal(pmns_squared_invariants, expected_squared_invariants, 1.e-10));
}


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cSMHdCKM_Kappa_running

#include "config.h"
#include "ew_input.hpp"
#include "lowe.h"

#include "cSMHdCKM_soft_parameters.hpp"

#include "matrix_tests.hpp"
#include "random_matrix.hpp"

#include <Eigen/LU>

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#ifdef ENABLE_RANDOM
#include <random>
#endif

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

void initialize_model_input(cSMHdCKM_input_parameters& input)
{
   input.LambdaIN = 0.2;
   input.Qin = 2.e16;
   input.QEWSB = Electroweak_constants::MTOP;
   input.sign_delta_mAsq = 1;

#ifdef ENABLE_RANDOM
   std::mt19937 generator;
   random_cue_matrix(input.UvInput, generator);
#endif
}

void initialize_model(cSMHdCKM_soft_parameters& model)
{
   model.set_scale(Electroweak_constants::MZ);
   model.set_loops(2);

   const double thetaw = ArcSin(Sqrt(0.231));
   const double cw = Cos(thetaw);
   const double sw = Sin(thetaw);
   const double v = 246.;

   model.set_g1(Sqrt(5. / 3.) * Sqrt(4. * Pi / (128. * cw * cw)));
   model.set_g2(Sqrt(4. * Pi / (128. * sw * sw)));
   model.set_g3(Sqrt(4. * Pi * 0.118));

   model.set_Yu(0, 0, 2.33e-3* Sqrt(2.) / v);
   model.set_Yu(1, 1, 0.677* Sqrt(2.) / v);
   model.set_Yu(2, 2, 172.7 * Sqrt(2.) / v);

   model.set_Yd(0, 0, 4.69e-3 * Sqrt(2.) / v);
   model.set_Yd(1, 1, 9.34e-2 * Sqrt(2.) / v);
   model.set_Yd(2, 2, 3.00 * Sqrt(2.) / v);

   model.set_Ye(0, 0, 4.8684727e-4 * Sqrt(2.) / v);
   model.set_Ye(1, 1, 1.0275138e-1 * Sqrt(2.) / v);
   model.set_Ye(2, 2, 1.7467 * Sqrt(2.) / v);

   model.set_Lambdax(0.25);

   model.set_mu2(Electroweak_constants::mu2SM);
   model.set_v(v);

   Eigen::Matrix<std::complex<double>, 3, 3> Kappa;
   Kappa << -3.32047e-15, -1.69446e-17, 1.69446e-17,
      -1.69446e-17, -4.02645e-15, -6.84641e-16,
      1.69446e-17, -6.84641e-16, -4.02645e-15;
   model.set_Kappa(Kappa);
}

BOOST_AUTO_TEST_CASE( test_Kappa_one_loop_running )
{
   cSMHdCKM_soft_parameters model;
   softsusy::QedQcd qedqcd;
   cSMHdCKM_input_parameters input;

   initialize_sm_input(qedqcd);
   initialize_model_input(input);
   initialize_model(model);

   model.set_loops(1);
   const double to_scale = 1000.;

   model.run_to(to_scale);

   // obtained using REAP-MPT version 1.9.3
   Eigen::Matrix<std::complex<double>, 3, 3> Kappa_expected;
   Kappa_expected << -3.54977e-15, -1.81147e-17, 1.81147e-17,
      -1.81147e-17, -4.3045e-15, -7.31917e-16,
      1.81147e-17, -7.31917e-16, -4.30448e-15;

   const auto inv_scale = Cbrt(Abs(Kappa_expected.determinant()));
   BOOST_CHECK(is_equal(model.get_Kappa() / inv_scale, Kappa_expected / inv_scale, 1.e-6));
}

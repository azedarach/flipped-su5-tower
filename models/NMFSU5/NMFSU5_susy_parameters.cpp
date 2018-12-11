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

#include "NMFSU5_susy_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME NMFSU5_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES(l) calc_susy_traces(l);

const int NMFSU5_susy_parameters::numberOfParameters;

NMFSU5_susy_parameters::NMFSU5_susy_parameters(const NMFSU5_input_parameters& input_)
   : input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

NMFSU5_susy_parameters::NMFSU5_susy_parameters(
   double scale_, int loops_, int thresholds_, const NMFSU5_input_parameters& input_,
   double g5_, double gX_, double lam1_, double lam2_, double lam3_, double lam3t_,
   double lam4_, double lam4t_, double lam5_, double lam5t_, double lam6_, double lam6t_,
   const std::complex<double>& lam7_, const std::complex<double>& lam8_,
   const std::complex<double>& eta1_, const std::complex<double>& eta2_,
   const std::complex<double>& eta3_, const Eigen::Matrix<std::complex<double>,3,3>& Y10_,
   const Eigen::Matrix<std::complex<double>,3,3>& Y10Pr_,
   const Eigen::Matrix<std::complex<double>,3,3>& Y5b_,
   const Eigen::Matrix<std::complex<double>,3,3>& Y5bPr_,
   const Eigen::Matrix<std::complex<double>,3,3>& Y1_,
   const Eigen::Matrix<std::complex<double>,3,3>& Y1Pr_
)
   : Beta_function()
   , g5(g5_), gX(gX_), lam1(lam1_), lam2(lam2_), lam3(lam3_), lam3t(lam3t_)
   , lam4(lam4_), lam4t(lam4t_), lam5(lam5_), lam5t(lam5t_)
   , lam6(lam6_), lam6t(lam6t_), lam7(lam7_), lam8(lam8_)
   , eta1(eta1_), eta2(eta2_), eta3(eta3_), Y10(Y10_), Y10Pr(Y10Pr_)
   , Y5b(Y5b_), Y5bPr(Y5bPr_), Y1(Y1_), Y1Pr(Y1Pr_)
   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd NMFSU5_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

NMFSU5_susy_parameters NMFSU5_susy_parameters::calc_beta(int loops) const
{
   double beta_g5 = 0.;
   double beta_gX = 0.;
   double beta_lam1 = 0.;
   double beta_lam2 = 0.;
   double beta_lam3 = 0.;
   double beta_lam3t = 0.;
   double beta_lam4 = 0.;
   double beta_lam4t = 0.;
   double beta_lam5 = 0.;
   double beta_lam5t = 0.;
   double beta_lam6 = 0.;
   double beta_lam6t = 0.;
   std::complex<double> beta_lam7 = 0.;
   std::complex<double> beta_lam8 = 0.;
   std::complex<double> beta_eta1 = 0.;
   std::complex<double> beta_eta2 = 0.;
   std::complex<double> beta_eta3 = 0.;
   Eigen::Matrix<std::complex<double>,3,3> beta_Y10 = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Y10Pr = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Y5b = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Y5bPr = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Y1 = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Y1Pr = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_g5 += calc_beta_g5_1_loop(TRACE_STRUCT);
      beta_gX += calc_beta_gX_1_loop(TRACE_STRUCT);
      beta_lam1 += calc_beta_lam1_1_loop(TRACE_STRUCT);
      beta_lam2 += calc_beta_lam2_1_loop(TRACE_STRUCT);
      beta_lam3 += calc_beta_lam3_1_loop(TRACE_STRUCT);
      beta_lam3t += calc_beta_lam3t_1_loop(TRACE_STRUCT);
      beta_lam4 += calc_beta_lam4_1_loop(TRACE_STRUCT);
      beta_lam4t += calc_beta_lam4t_1_loop(TRACE_STRUCT);
      beta_lam5 += calc_beta_lam5_1_loop(TRACE_STRUCT);
      beta_lam5t += calc_beta_lam5t_1_loop(TRACE_STRUCT);
      beta_lam6 += calc_beta_lam6_1_loop(TRACE_STRUCT);
      beta_lam6t += calc_beta_lam6t_1_loop(TRACE_STRUCT);
      beta_lam7 += calc_beta_lam7_1_loop(TRACE_STRUCT);
      beta_lam8 += calc_beta_lam8_1_loop(TRACE_STRUCT);
      beta_eta1 += calc_beta_eta1_1_loop(TRACE_STRUCT);
      beta_eta2 += calc_beta_eta2_1_loop(TRACE_STRUCT);
      beta_eta3 += calc_beta_eta3_1_loop(TRACE_STRUCT);
      beta_Y10 += calc_beta_Y10_1_loop(TRACE_STRUCT);
      beta_Y10Pr += calc_beta_Y10Pr_1_loop(TRACE_STRUCT);
      beta_Y5b += calc_beta_Y5b_1_loop(TRACE_STRUCT);
      beta_Y5bPr += calc_beta_Y5bPr_1_loop(TRACE_STRUCT);
      beta_Y1 += calc_beta_Y1_1_loop(TRACE_STRUCT);
      beta_Y1Pr += calc_beta_Y1Pr_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_g5 += calc_beta_g5_2_loop(TRACE_STRUCT);
         beta_gX += calc_beta_gX_2_loop(TRACE_STRUCT);
         beta_lam1 += calc_beta_lam1_2_loop(TRACE_STRUCT);
         beta_lam2 += calc_beta_lam2_2_loop(TRACE_STRUCT);
         beta_lam3 += calc_beta_lam3_2_loop(TRACE_STRUCT);
         beta_lam3t += calc_beta_lam3t_2_loop(TRACE_STRUCT);
         beta_lam4 += calc_beta_lam4_2_loop(TRACE_STRUCT);
         beta_lam4t += calc_beta_lam4t_2_loop(TRACE_STRUCT);
         beta_lam5 += calc_beta_lam5_2_loop(TRACE_STRUCT);
         beta_lam5t += calc_beta_lam5t_2_loop(TRACE_STRUCT);
         beta_lam6 += calc_beta_lam6_2_loop(TRACE_STRUCT);
         beta_lam6t += calc_beta_lam6t_2_loop(TRACE_STRUCT);
         beta_lam7 += calc_beta_lam7_2_loop(TRACE_STRUCT);
         beta_lam8 += calc_beta_lam8_2_loop(TRACE_STRUCT);
         beta_eta1 += calc_beta_eta1_2_loop(TRACE_STRUCT);
         beta_eta2 += calc_beta_eta2_2_loop(TRACE_STRUCT);
         beta_eta3 += calc_beta_eta3_2_loop(TRACE_STRUCT);
         beta_Y10 += calc_beta_Y10_2_loop(TRACE_STRUCT);
         beta_Y10Pr += calc_beta_Y10Pr_2_loop(TRACE_STRUCT);
         beta_Y5b += calc_beta_Y5b_2_loop(TRACE_STRUCT);
         beta_Y5bPr += calc_beta_Y5bPr_2_loop(TRACE_STRUCT);
         beta_Y1 += calc_beta_Y1_2_loop(TRACE_STRUCT);
         beta_Y1Pr += calc_beta_Y1Pr_2_loop(TRACE_STRUCT);

         if (loops > 3) {
#ifdef ENABLE_THREADS
            {
               auto fut_g5 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_g5_3_loop(TRACE_STRUCT); });
               auto fut_gX = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_gX_3_loop(TRACE_STRUCT); });
               auto fut_lam1 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam1_3_loop(TRACE_STRUCT); });
               auto fut_lam2 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam2_3_loop(TRACE_STRUCT); });
               auto fut_lam3 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam3_3_loop(TRACE_STRUCT); });
               auto fut_lam3t = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam3t_3_loop(TRACE_STRUCT); });
               auto fut_lam4 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam4_3_loop(TRACE_STRUCT); });
               auto fut_lam4t = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam4t_3_loop(TRACE_STRUCT); });
               auto fut_lam5 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam5_3_loop(TRACE_STRUCT); });
               auto fut_lam5t = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam5t_3_loop(TRACE_STRUCT); });
               auto fut_lam6 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam6_3_loop(TRACE_STRUCT); });
               auto fut_lam6t = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam6t_3_loop(TRACE_STRUCT); });
               auto fut_lam7 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam7_3_loop(TRACE_STRUCT); });
               auto fut_lam8 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_lam8_3_loop(TRACE_STRUCT); });
               auto fut_eta1 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_eta1_3_loop(TRACE_STRUCT); });
               auto fut_eta2 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_eta2_3_loop(TRACE_STRUCT); });
               auto fut_eta3 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_eta3_3_loop(TRACE_STRUCT); });
               auto fut_Y10 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y10_3_loop(TRACE_STRUCT); });
               auto fut_Y10Pr = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y10Pr_3_loop(TRACE_STRUCT); });
               auto fut_Y5b = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y5b_3_loop(TRACE_STRUCT); });
               auto fut_Y5bPr = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y5bPr_3_loop(TRACE_STRUCT); });
               auto fut_Y1 = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y1_3_loop(TRACE_STRUCT); });
               auto fut_Y1Pr = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_Y1Pr_3_loop(TRACE_STRUCT); });

               beta_g5 += fut_g5.get();
               beta_gX += fut_gX.get();
               beta_lam1 += fut_lam1.get();
               beta_lam2 += fut_lam2.get();
               beta_lam3 += fut_lam3.get();
               beta_lam3t += fut_lam3t.get();
               beta_lam4 += fut_lam4.get();
               beta_lam4t += fut_lam4t.get();
               beta_lam5 += fut_lam5.get();
               beta_lam5t += fut_lam5t.get();
               beta_lam6 += fut_lam6.get();
               beta_lam6t += fut_lam6t.get();
               beta_lam7 += fut_lam7.get();
               beta_lam8 += fut_lam8.get();
               beta_eta1 += fut_eta1.get();
               beta_eta2 += fut_eta2.get();
               beta_eta3 += fut_eta3.get();
               beta_Y10 += fut_Y10.get();
               beta_Y10Pr += fut_Y10Pr.get();
               beta_Y5b += fut_Y5b.get();
               beta_Y5bPr += fut_Y5bPr.get();
               beta_Y1 += fut_Y1.get();
               beta_Y1Pr += fut_Y1Pr.get();
            }
#else
            beta_g5 += calc_beta_g5_3_loop(TRACE_STRUCT);
            beta_gX += calc_beta_gX_3_loop(TRACE_STRUCT);
            beta_lam1 += calc_beta_lam1_3_loop(TRACE_STRUCT);
            beta_lam2 += calc_beta_lam2_3_loop(TRACE_STRUCT);
            beta_lam3 += calc_beta_lam3_3_loop(TRACE_STRUCT);
            beta_lam3t += calc_beta_lam3t_3_loop(TRACE_STRUCT);
            beta_lam4 += calc_beta_lam4_3_loop(TRACE_STRUCT);
            beta_lam4t += calc_beta_lam4t_3_loop(TRACE_STRUCT);
            beta_lam5 += calc_beta_lam5_3_loop(TRACE_STRUCT);
            beta_lam5t += calc_beta_lam5t_3_loop(TRACE_STRUCT);
            beta_lam6 += calc_beta_lam6_3_loop(TRACE_STRUCT);
            beta_lam6t += calc_beta_lam6t_3_loop(TRACE_STRUCT);
            beta_lam7 += calc_beta_lam7_3_loop(TRACE_STRUCT);
            beta_lam8 += calc_beta_lam8_3_loop(TRACE_STRUCT);
            beta_eta1 += calc_beta_eta1_3_loop(TRACE_STRUCT);
            beta_eta2 += calc_beta_eta2_3_loop(TRACE_STRUCT);
            beta_eta3 += calc_beta_eta3_3_loop(TRACE_STRUCT);
            beta_Y10 += calc_beta_Y10_3_loop(TRACE_STRUCT);
            beta_Y10Pr += calc_beta_Y10Pr_3_loop(TRACE_STRUCT);
            beta_Y5b += calc_beta_Y5b_3_loop(TRACE_STRUCT);
            beta_Y5bPr += calc_beta_Y5bPr_3_loop(TRACE_STRUCT);
            beta_Y1 += calc_beta_Y1_3_loop(TRACE_STRUCT);
            beta_Y1Pr += calc_beta_Y1Pr_3_loop(TRACE_STRUCT);
#endif

            if (loops > 3) {
               beta_g5 += calc_beta_g5_4_loop(TRACE_STRUCT);
               beta_gX += calc_beta_gX_4_loop(TRACE_STRUCT);
               beta_lam1 += calc_beta_lam1_4_loop(TRACE_STRUCT);
               beta_lam2 += calc_beta_lam2_4_loop(TRACE_STRUCT);
               beta_lam3 += calc_beta_lam3_4_loop(TRACE_STRUCT);
               beta_lam3t += calc_beta_lam3t_4_loop(TRACE_STRUCT);
               beta_lam4 += calc_beta_lam4_4_loop(TRACE_STRUCT);
               beta_lam4t += calc_beta_lam4t_4_loop(TRACE_STRUCT);
               beta_lam5 += calc_beta_lam5_4_loop(TRACE_STRUCT);
               beta_lam5t += calc_beta_lam5t_4_loop(TRACE_STRUCT);
               beta_lam6 += calc_beta_lam6_4_loop(TRACE_STRUCT);
               beta_lam6t += calc_beta_lam6t_4_loop(TRACE_STRUCT);
               beta_lam7 += calc_beta_lam7_4_loop(TRACE_STRUCT);
               beta_lam8 += calc_beta_lam8_4_loop(TRACE_STRUCT);
               beta_eta1 += calc_beta_eta1_4_loop(TRACE_STRUCT);
               beta_eta2 += calc_beta_eta2_4_loop(TRACE_STRUCT);
               beta_eta3 += calc_beta_eta3_4_loop(TRACE_STRUCT);
               beta_Y10 += calc_beta_Y10_4_loop(TRACE_STRUCT);
               beta_Y10Pr += calc_beta_Y10Pr_4_loop(TRACE_STRUCT);
               beta_Y5b += calc_beta_Y5b_4_loop(TRACE_STRUCT);
               beta_Y5bPr += calc_beta_Y5bPr_4_loop(TRACE_STRUCT);
               beta_Y1 += calc_beta_Y1_4_loop(TRACE_STRUCT);
               beta_Y1Pr += calc_beta_Y1Pr_4_loop(TRACE_STRUCT);
            }
         }
      }
   }

   return NMFSU5_susy_parameters(get_scale(), loops, get_thresholds(), input,
                                 beta_g5, beta_gX, beta_lam1, beta_lam2, beta_lam3,
                                 beta_lam3t, beta_lam4, beta_lam4t, beta_lam5,
                                 beta_lam5t, beta_lam6, beta_lam6t,
                                 beta_lam7, beta_lam8, beta_eta1, beta_eta2,
                                 beta_eta3, beta_Y10, beta_Y10Pr, beta_Y5b,
                                 beta_Y5bPr, beta_Y1, beta_Y1Pr);

}

NMFSU5_susy_parameters NMFSU5_susy_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void NMFSU5_susy_parameters::clear()
{
   reset();
   g5 = 0.;
   gX = 0.;
   lam1 = 0.;
   lam2 = 0.;
   lam3 = 0.;
   lam3t = 0.;
   lam4 = 0.;
   lam4t = 0.;
   lam5 = 0.;
   lam5t = 0.;
   lam6 = 0.;
   lam6t = 0.;
   lam7 = 0.;
   lam8 = 0.;
   eta1 = 0.;
   eta2 = 0.;
   eta3 = 0.;
   Y10 = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Y10Pr = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Y5b = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Y5bPr = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Y1 = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Y1Pr = Eigen::Matrix<std::complex<double>,3,3>::Zero();
}

Eigen::ArrayXd NMFSU5_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g5;
   pars(1) = gX;
   pars(2) = lam1;
   pars(3) = lam2;
   pars(4) = lam3;
   pars(5) = lam3t;
   pars(6) = lam4;
   pars(7) = lam4t;
   pars(8) = lam5;
   pars(9) = lam5t;
   pars(10) = lam6;
   pars(11) = lam6t;
   pars(12) = Re(lam7);
   pars(13) = Im(lam7);
   pars(14) = Re(lam8);
   pars(15) = Im(lam8);
   pars(16) = Re(eta1);
   pars(17) = Im(eta1);
   pars(18) = Re(eta2);
   pars(19) = Im(eta2);
   pars(20) = Re(eta3);
   pars(21) = Im(eta3);
   pars(22) = Re(Y10(0,0));
   pars(23) = Im(Y10(0,0));
   pars(24) = Re(Y10(0,1));
   pars(25) = Im(Y10(0,1));
   pars(26) = Re(Y10(0,2));
   pars(27) = Im(Y10(0,2));
   pars(28) = Re(Y10(1,0));
   pars(29) = Im(Y10(1,0));
   pars(30) = Re(Y10(1,1));
   pars(31) = Im(Y10(1,1));
   pars(32) = Re(Y10(1,2));
   pars(33) = Im(Y10(1,2));
   pars(34) = Re(Y10(2,0));
   pars(35) = Im(Y10(2,0));
   pars(36) = Re(Y10(2,1));
   pars(37) = Im(Y10(2,1));
   pars(38) = Re(Y10(2,2));
   pars(39) = Im(Y10(2,2));
   pars(40) = Re(Y10Pr(0,0));
   pars(41) = Im(Y10Pr(0,0));
   pars(42) = Re(Y10Pr(0,1));
   pars(43) = Im(Y10Pr(0,1));
   pars(44) = Re(Y10Pr(0,2));
   pars(45) = Im(Y10Pr(0,2));
   pars(46) = Re(Y10Pr(1,0));
   pars(47) = Im(Y10Pr(1,0));
   pars(48) = Re(Y10Pr(1,1));
   pars(49) = Im(Y10Pr(1,1));
   pars(50) = Re(Y10Pr(1,2));
   pars(51) = Im(Y10Pr(1,2));
   pars(52) = Re(Y10Pr(2,0));
   pars(53) = Im(Y10Pr(2,0));
   pars(54) = Re(Y10Pr(2,1));
   pars(55) = Im(Y10Pr(2,1));
   pars(56) = Re(Y10Pr(2,2));
   pars(57) = Im(Y10Pr(2,2));
   pars(58) = Re(Y5b(0,0));
   pars(59) = Im(Y5b(0,0));
   pars(60) = Re(Y5b(0,1));
   pars(61) = Im(Y5b(0,1));
   pars(62) = Re(Y5b(0,2));
   pars(63) = Im(Y5b(0,2));
   pars(64) = Re(Y5b(1,0));
   pars(65) = Im(Y5b(1,0));
   pars(66) = Re(Y5b(1,1));
   pars(67) = Im(Y5b(1,1));
   pars(68) = Re(Y5b(1,2));
   pars(69) = Im(Y5b(1,2));
   pars(70) = Re(Y5b(2,0));
   pars(71) = Im(Y5b(2,0));
   pars(72) = Re(Y5b(2,1));
   pars(73) = Im(Y5b(2,1));
   pars(74) = Re(Y5b(2,2));
   pars(75) = Im(Y5b(2,2));
   pars(76) = Re(Y5bPr(0,0));
   pars(77) = Im(Y5bPr(0,0));
   pars(78) = Re(Y5bPr(0,1));
   pars(79) = Im(Y5bPr(0,1));
   pars(80) = Re(Y5bPr(0,2));
   pars(81) = Im(Y5bPr(0,2));
   pars(82) = Re(Y5bPr(1,0));
   pars(83) = Im(Y5bPr(1,0));
   pars(84) = Re(Y5bPr(1,1));
   pars(85) = Im(Y5bPr(1,1));
   pars(86) = Re(Y5bPr(1,2));
   pars(87) = Im(Y5bPr(1,2));
   pars(88) = Re(Y5bPr(2,0));
   pars(89) = Im(Y5bPr(2,0));
   pars(90) = Re(Y5bPr(2,1));
   pars(91) = Im(Y5bPr(2,1));
   pars(92) = Re(Y5bPr(2,2));
   pars(93) = Im(Y5bPr(2,2));
   pars(94) = Re(Y1(0,0));
   pars(95) = Im(Y1(0,0));
   pars(96) = Re(Y1(0,1));
   pars(97) = Im(Y1(0,1));
   pars(98) = Re(Y1(0,2));
   pars(99) = Im(Y1(0,2));
   pars(100) = Re(Y1(1,0));
   pars(101) = Im(Y1(1,0));
   pars(102) = Re(Y1(1,1));
   pars(103) = Im(Y1(1,1));
   pars(104) = Re(Y1(1,2));
   pars(105) = Im(Y1(1,2));
   pars(106) = Re(Y1(2,0));
   pars(107) = Im(Y1(2,0));
   pars(108) = Re(Y1(2,1));
   pars(109) = Im(Y1(2,1));
   pars(110) = Re(Y1(2,2));
   pars(111) = Im(Y1(2,2));
   pars(112) = Re(Y1Pr(0,0));
   pars(113) = Im(Y1Pr(0,0));
   pars(114) = Re(Y1Pr(0,1));
   pars(115) = Im(Y1Pr(0,1));
   pars(116) = Re(Y1Pr(0,2));
   pars(117) = Im(Y1Pr(0,2));
   pars(118) = Re(Y1Pr(1,0));
   pars(119) = Im(Y1Pr(1,0));
   pars(120) = Re(Y1Pr(1,1));
   pars(121) = Im(Y1Pr(1,1));
   pars(122) = Re(Y1Pr(1,2));
   pars(123) = Im(Y1Pr(1,2));
   pars(124) = Re(Y1Pr(2,0));
   pars(125) = Im(Y1Pr(2,0));
   pars(126) = Re(Y1Pr(2,1));
   pars(127) = Im(Y1Pr(2,1));
   pars(128) = Re(Y1Pr(2,2));
   pars(129) = Im(Y1Pr(2,2));

   return pars;

}

void NMFSU5_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "g5 = " << g5 << '\n';
   ostr << "gX = " << gX << '\n';
   ostr << "lam1 = " << lam1 << '\n';
   ostr << "lam2 = " << lam2 << '\n';
   ostr << "lam3 = " << lam3 << '\n';
   ostr << "lam3t = " << lam3t << '\n';
   ostr << "lam4 = " << lam4 << '\n';
   ostr << "lam4t = " << lam4t << '\n';
   ostr << "lam5 = " << lam5 << '\n';
   ostr << "lam5t = " << lam5t << '\n';
   ostr << "lam6 = " << lam6 << '\n';
   ostr << "lam6t = " << lam6t << '\n';
   ostr << "lam7 = " << lam7 << '\n';
   ostr << "lam8 = " << lam8 << '\n';
   ostr << "eta1 = " << eta1 << '\n';
   ostr << "eta2 = " << eta2 << '\n';
   ostr << "eta3 = " << eta3 << '\n';
   ostr << "Y10 = " << Y10 << '\n';
   ostr << "Y10Pr = " << Y10Pr << '\n';
   ostr << "Y5b = " << Y5b << '\n';
   ostr << "Y5bPr = " << Y5bPr << '\n';
   ostr << "Y1 = " << Y1 << '\n';
   ostr << "Y1Pr = " << Y1Pr << '\n';
}

void NMFSU5_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   g5 = pars(0);
   gX = pars(1);
   lam1 = pars(2);
   lam2 = pars(3);
   lam3 = pars(4);
   lam3t = pars(5);
   lam4 = pars(6);
   lam4t = pars(7);
   lam5 = pars(8);
   lam5t = pars(9);
   lam6 = pars(10);
   lam6t = pars(11);
   lam7 = std::complex<double>(pars(12), pars(13));
   lam8 = std::complex<double>(pars(14), pars(15));
   eta1 = std::complex<double>(pars(16), pars(17));
   eta2 = std::complex<double>(pars(18), pars(19));
   eta3 = std::complex<double>(pars(20), pars(21));
   Y10(0,0) = std::complex<double>(pars(22), pars(23));
   Y10(0,1) = std::complex<double>(pars(24), pars(25));
   Y10(0,2) = std::complex<double>(pars(26), pars(27));
   Y10(1,0) = std::complex<double>(pars(28), pars(29));
   Y10(1,1) = std::complex<double>(pars(30), pars(31));
   Y10(1,2) = std::complex<double>(pars(32), pars(33));
   Y10(2,0) = std::complex<double>(pars(34), pars(35));
   Y10(2,1) = std::complex<double>(pars(36), pars(37));
   Y10(2,2) = std::complex<double>(pars(38), pars(39));
   Y10Pr(0,0) = std::complex<double>(pars(40), pars(41));
   Y10Pr(0,1) = std::complex<double>(pars(42), pars(43));
   Y10Pr(0,2) = std::complex<double>(pars(44), pars(45));
   Y10Pr(1,0) = std::complex<double>(pars(46), pars(47));
   Y10Pr(1,1) = std::complex<double>(pars(48), pars(49));
   Y10Pr(1,2) = std::complex<double>(pars(50), pars(51));
   Y10Pr(2,0) = std::complex<double>(pars(52), pars(53));
   Y10Pr(2,1) = std::complex<double>(pars(54), pars(55));
   Y10Pr(2,2) = std::complex<double>(pars(56), pars(57));
   Y5b(0,0) = std::complex<double>(pars(58), pars(59));
   Y5b(0,1) = std::complex<double>(pars(60), pars(61));
   Y5b(0,2) = std::complex<double>(pars(62), pars(63));
   Y5b(1,0) = std::complex<double>(pars(64), pars(65));
   Y5b(1,1) = std::complex<double>(pars(66), pars(67));
   Y5b(1,2) = std::complex<double>(pars(68), pars(69));
   Y5b(2,0) = std::complex<double>(pars(70), pars(71));
   Y5b(2,1) = std::complex<double>(pars(72), pars(73));
   Y5b(2,2) = std::complex<double>(pars(74), pars(75));
   Y5bPr(0,0) = std::complex<double>(pars(76), pars(77));
   Y5bPr(0,1) = std::complex<double>(pars(78), pars(79));
   Y5bPr(0,2) = std::complex<double>(pars(80), pars(81));
   Y5bPr(1,0) = std::complex<double>(pars(82), pars(83));
   Y5bPr(1,1) = std::complex<double>(pars(84), pars(85));
   Y5bPr(1,2) = std::complex<double>(pars(86), pars(87));
   Y5bPr(2,0) = std::complex<double>(pars(88), pars(89));
   Y5bPr(2,1) = std::complex<double>(pars(90), pars(91));
   Y5bPr(2,2) = std::complex<double>(pars(92), pars(93));
   Y1(0,0) = std::complex<double>(pars(94), pars(95));
   Y1(0,1) = std::complex<double>(pars(96), pars(97));
   Y1(0,2) = std::complex<double>(pars(98), pars(99));
   Y1(1,0) = std::complex<double>(pars(100), pars(101));
   Y1(1,1) = std::complex<double>(pars(102), pars(103));
   Y1(1,2) = std::complex<double>(pars(104), pars(105));
   Y1(2,0) = std::complex<double>(pars(106), pars(107));
   Y1(2,1) = std::complex<double>(pars(108), pars(109));
   Y1(2,2) = std::complex<double>(pars(110), pars(111));
   Y1Pr(0,0) = std::complex<double>(pars(112), pars(113));
   Y1Pr(0,1) = std::complex<double>(pars(114), pars(115));
   Y1Pr(0,2) = std::complex<double>(pars(116), pars(117));
   Y1Pr(1,0) = std::complex<double>(pars(118), pars(119));
   Y1Pr(1,1) = std::complex<double>(pars(120), pars(121));
   Y1Pr(1,2) = std::complex<double>(pars(122), pars(123));
   Y1Pr(2,0) = std::complex<double>(pars(124), pars(125));
   Y1Pr(2,1) = std::complex<double>(pars(126), pars(127));
   Y1Pr(2,2) = std::complex<double>(pars(128), pars(129));

}

const NMFSU5_input_parameters& NMFSU5_susy_parameters::get_input() const
{
   return input;
}

NMFSU5_input_parameters& NMFSU5_susy_parameters::get_input()
{
   return input;
}

void NMFSU5_susy_parameters::set_input_parameters(const NMFSU5_input_parameters& input_)
{
   input = input_;
}

NMFSU5_susy_parameters::Susy_traces NMFSU5_susy_parameters::calc_susy_traces(int loops) const
{
   Susy_traces susy_traces;

   if (loops > 0) {

   }

   if (loops > 1) {

   }

   if (loops > 2) {

   }

   return susy_traces;
}

std::ostream& operator<<(std::ostream& ostr, const NMFSU5_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

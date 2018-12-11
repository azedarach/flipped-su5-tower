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

#include "NMFSU5_soft_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <Eigen/LU>

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int NMFSU5_soft_parameters::numberOfParameters;

NMFSU5_soft_parameters::NMFSU5_soft_parameters(const NMFSU5_input_parameters& input_)
   : NMFSU5_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

NMFSU5_soft_parameters::NMFSU5_soft_parameters(
   const NMFSU5_susy_parameters& susy_model, double m10sq_, double m5sq_, double m5Prsq_,
   const std::complex<double>& m12sq_, const std::complex<double>& mu_,
   const std::complex<double>& muPr_, double v_, double vPr_, double VG_)
   : NMFSU5_susy_parameters(susy_model)
   , m10sq(m10sq_), m5sq(m5sq_), m5Prsq(m5Prsq_), m12sq(m12sq_)
   , mu(mu_), muPr(muPr_), v(v_), vPr(vPr_), VG(VG_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd NMFSU5_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

NMFSU5_soft_parameters NMFSU5_soft_parameters::calc_beta(int loops) const
{
   double beta_m10sq = 0.;
   double beta_m5sq = 0.;
   double beta_m5Prsq = 0.;
   std::complex<double> beta_m12sq = 0.;
   std::complex<double> beta_mu = 0.;
   std::complex<double> beta_muPr = 0.;
   double beta_v = 0.;
   double beta_vPr = 0.;
   double beta_VG = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_m10sq += calc_beta_m10sq_1_loop(TRACE_STRUCT);
      beta_m5sq += calc_beta_m5sq_1_loop(TRACE_STRUCT);
      beta_m5Prsq += calc_beta_m5Prsq_1_loop(TRACE_STRUCT);
      beta_m12sq += calc_beta_m12sq_1_loop(TRACE_STRUCT);
      beta_mu += calc_beta_mu_1_loop(TRACE_STRUCT);
      beta_muPr += calc_beta_muPr_1_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_1_loop(TRACE_STRUCT);
      beta_vPr += calc_beta_vPr_1_loop(TRACE_STRUCT);
      beta_VG += calc_beta_VG_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_m10sq += calc_beta_m10sq_2_loop(TRACE_STRUCT);
         beta_m5sq += calc_beta_m5sq_2_loop(TRACE_STRUCT);
         beta_m5Prsq += calc_beta_m5Prsq_2_loop(TRACE_STRUCT);
         beta_m12sq += calc_beta_m12sq_2_loop(TRACE_STRUCT);
         beta_mu += calc_beta_mu_2_loop(TRACE_STRUCT);
         beta_muPr += calc_beta_muPr_2_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_2_loop(TRACE_STRUCT);
         beta_vPr += calc_beta_vPr_2_loop(TRACE_STRUCT);
         beta_VG += calc_beta_VG_2_loop(TRACE_STRUCT);

         if (loops > 2) {
#ifdef ENABLE_THREADS
            {
               auto fut_m10sq = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_m10sq_3_loop(TRACE_STRUCT); });
               auto fut_m5sq = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_m5sq_3_loop(TRACE_STRUCT); });
               auto fut_m5Prsq = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_m5Prsq_3_loop(TRACE_STRUCT); });
               auto fut_m12sq = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_m12sq_3_loop(TRACE_STRUCT); });
               auto fut_mu = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_mu_3_loop(TRACE_STRUCT); });
               auto fut_muPr = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_muPr_3_loop(TRACE_STRUCT); });
               auto fut_v = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_v_3_loop(TRACE_STRUCT); });
               auto fut_vPr = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_vPr_3_loop(TRACE_STRUCT); });
               auto fut_VG = global_thread_pool().run_packaged_task(
                  [this, &TRACE_STRUCT]() { return calc_beta_VG_3_loop(TRACE_STRUCT); });

               beta_m10sq += fut_m10sq.get();
               beta_m5sq += fut_m5sq.get();
               beta_m5Prsq += fut_m5Prsq.get();
               beta_m12sq += fut_m12sq.get();
               beta_mu += fut_mu.get();
               beta_muPr += fut_muPr.get();
               beta_v += fut_v.get();
               beta_vPr += fut_vPr.get();
               beta_VG += fut_VG.get();
            }
#else
            beta_m10sq += calc_beta_m10sq_3_loop(TRACE_STRUCT);
            beta_m5sq += calc_beta_m5sq_3_loop(TRACE_STRUCT);
            beta_m5Prsq += calc_beta_m5Prsq_3_loop(TRACE_STRUCT);
            beta_m12sq += calc_beta_m12sq_3_loop(TRACE_STRUCT);
            beta_mu += calc_beta_mu_3_loop(TRACE_STRUCT);
            beta_muPr += calc_beta_muPr_3_loop(TRACE_STRUCT);
            beta_v += calc_beta_v_3_loop(TRACE_STRUCT);
            beta_vPr += calc_beta_vPr_3_loop(TRACE_STRUCT);
            beta_VG += calc_beta_VG_3_loop(TRACE_STRUCT);
#endif

            if (loops > 3) {

            }
         }
      }
   }

   const NMFSU5_susy_parameters susy_betas(NMFSU5_susy_parameters::calc_beta(loops));

   return NMFSU5_soft_parameters(susy_betas, beta_m10sq, beta_m5sq, beta_m5Prsq,
                                 beta_m12sq, beta_mu, beta_muPr, beta_v, beta_vPr, beta_VG);
}

NMFSU5_soft_parameters NMFSU5_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void NMFSU5_soft_parameters::clear()
{
   NMFSU5_susy_parameters::clear();

   m10sq = 0.;
   m5sq = 0.;
   m5Prsq = 0.;
   m12sq = 0.;
   mu = 0.;
   muPr = 0.;
   v = 0.;
   vPr = 0.;
   VG = 0.;
}

Eigen::ArrayXd NMFSU5_soft_parameters::get() const
{
   Eigen::ArrayXd pars(NMFSU5_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(130) = m10sq;
   pars(131) = m5sq;
   pars(132) = m5Prsq;
   pars(133) = Re(m12sq);
   pars(134) = Im(m12sq);
   pars(135) = Re(mu);
   pars(136) = Im(mu);
   pars(137) = Re(muPr);
   pars(138) = Im(muPr);
   pars(139) = v;
   pars(140) = vPr;
   pars(141) = VG;

   return pars;
}

void NMFSU5_soft_parameters::print(std::ostream& ostr) const
{
   NMFSU5_susy_parameters::print(ostr);
   ostr<< "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "m10sq = " << m10sq << '\n';
   ostr << "m5sq = " << m5sq << '\n';
   ostr << "m5Prsq = " << m5Prsq << '\n';
   ostr << "m12sq = " << m12sq << '\n';
   ostr << "mu = " << mu << '\n';
   ostr << "muPr = " << muPr << '\n';
   ostr << "v = " << v << '\n';
   ostr << "vPr = " << vPr << '\n';
   ostr << "VG = " << VG << '\n';
}

void NMFSU5_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   NMFSU5_susy_parameters::set(pars);

   m10sq = pars(130);
   m5sq = pars(131);
   m5Prsq = pars(132);
   m12sq = std::complex<double>(pars(133), pars(134));
   mu = std::complex<double>(pars(135), pars(136));
   muPr = std::complex<double>(pars(137), pars(138));
   v = pars(139);
   vPr = pars(140);
   VG = pars(141);
}

NMFSU5_soft_parameters::Soft_traces NMFSU5_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {

   }

   if (loops > 1) {

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const NMFSU5_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

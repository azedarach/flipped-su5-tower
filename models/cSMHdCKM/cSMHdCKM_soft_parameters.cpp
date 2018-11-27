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

// File generated at Mon 5 Nov 2018 12:47:13

#include "cSMHdCKM_soft_parameters.hpp"
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

const int cSMHdCKM_soft_parameters::numberOfParameters;

cSMHdCKM_soft_parameters::cSMHdCKM_soft_parameters(const cSMHdCKM_input_parameters& input_)
   : cSMHdCKM_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

cSMHdCKM_soft_parameters::cSMHdCKM_soft_parameters(
   const cSMHdCKM_susy_parameters& susy_model
   , double mu2_, double v_, double Kappa_inv_scale_,
   const Eigen::Matrix<std::complex<double>,3,3>& KappaND_
)
   : cSMHdCKM_susy_parameters(susy_model)
   , mu2(mu2_), v(v_), Kappa_inv_scale(Abs(Kappa_inv_scale_))
   , KappaND(KappaND_)
{
   set_number_of_parameters(numberOfParameters);
}

cSMHdCKM_soft_parameters::cSMHdCKM_soft_parameters(
   const cSMHdCKM_susy_parameters& susy_model
   , double mu2_, double v_, const Eigen::Matrix<std::complex<double>,3,3>& Kappa_
)
   : cSMHdCKM_susy_parameters(susy_model)
   , mu2(mu2_), v(v_)
{
   set_number_of_parameters(numberOfParameters);
   calculate_KappaND(Kappa_);
}

void cSMHdCKM_soft_parameters::calculate_KappaND(const Eigen::Matrix<std::complex<double>,3,3>& Kappa_)
{
   Kappa_inv_scale = Cbrt(Abs(Kappa_.determinant()));
   if (is_zero(Kappa_inv_scale)) {
      Kappa_inv_scale = Abs(Kappa_.mean());
   }
   if (!is_zero(Kappa_inv_scale)) {
      KappaND = Kappa_ / Kappa_inv_scale;
   } else {
      Kappa_inv_scale = 1.;
      KappaND = Kappa_;
   }
}

void cSMHdCKM_soft_parameters::set_Kappa(const Eigen::Matrix<std::complex<double>,3,3>& Kappa_)
{
   calculate_KappaND(Kappa_);
}

void cSMHdCKM_soft_parameters::set_Kappa(int i, int k, const std::complex<double>& value)
{
   Eigen::Matrix<std::complex<double>,3,3> Kappa(get_Kappa());
   Kappa(i,k) = value;
   calculate_KappaND(Kappa);
}

Eigen::ArrayXd cSMHdCKM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

cSMHdCKM_soft_parameters cSMHdCKM_soft_parameters::calc_beta(int loops) const
{
   double beta_mu2 = 0.;
   double beta_v = 0.;
   Eigen::Matrix<std::complex<double>,3,3> beta_KappaND
      = Eigen::Matrix<std::complex<double>,3,3>::Zero();

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_1_loop(TRACE_STRUCT);
      beta_KappaND += calc_beta_KappaND_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_2_loop(TRACE_STRUCT);
         beta_KappaND += calc_beta_KappaND_2_loop(TRACE_STRUCT);

         if (loops > 2) {
         #ifdef ENABLE_THREADS
            {
               auto fut_mu2 = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_mu2_3_loop(TRACE_STRUCT); });

               beta_mu2 += fut_mu2.get();

            }
         #else
            beta_mu2 += calc_beta_mu2_3_loop(TRACE_STRUCT);
         #endif

            if (loops > 3) {

            }
         }
      }
   }

   const cSMHdCKM_susy_parameters susy_betas(cSMHdCKM_susy_parameters::calc_beta(loops));

   return cSMHdCKM_soft_parameters(susy_betas, beta_mu2, beta_v, Kappa_inv_scale, beta_KappaND);
}

cSMHdCKM_soft_parameters cSMHdCKM_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void cSMHdCKM_soft_parameters::clear()
{
   cSMHdCKM_susy_parameters::clear();

   mu2 = 0.;
   v = 0.;
   Kappa_inv_scale = 1.;
   KappaND = Eigen::Matrix<std::complex<double>,3,3>::Zero();
}

Eigen::ArrayXd cSMHdCKM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(cSMHdCKM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(58) = mu2;
   pars(59) = v;
   pars(60) = Re(KappaND(0,0));
   pars(61) = Im(KappaND(0,0));
   pars(62) = Re(KappaND(0,1));
   pars(63) = Im(KappaND(0,1));
   pars(64) = Re(KappaND(0,2));
   pars(65) = Im(KappaND(0,2));
   pars(66) = Re(KappaND(1,0));
   pars(67) = Im(KappaND(1,0));
   pars(68) = Re(KappaND(1,1));
   pars(69) = Im(KappaND(1,1));
   pars(70) = Re(KappaND(1,2));
   pars(71) = Im(KappaND(1,2));
   pars(72) = Re(KappaND(2,0));
   pars(73) = Im(KappaND(2,0));
   pars(74) = Re(KappaND(2,1));
   pars(75) = Im(KappaND(2,1));
   pars(76) = Re(KappaND(2,2));
   pars(77) = Im(KappaND(2,2));

   return pars;
}

void cSMHdCKM_soft_parameters::print(std::ostream& ostr) const
{
   cSMHdCKM_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';
   ostr << "Kappa = " << Kappa_inv_scale * KappaND << '\n';
}

void cSMHdCKM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   cSMHdCKM_susy_parameters::set(pars);

   mu2 = pars(58);
   v = pars(59);
   KappaND(0,0) = std::complex<double>(pars(60), pars(61));
   KappaND(0,1) = std::complex<double>(pars(62), pars(63));
   KappaND(0,2) = std::complex<double>(pars(64), pars(65));
   KappaND(1,0) = std::complex<double>(pars(66), pars(67));
   KappaND(1,1) = std::complex<double>(pars(68), pars(69));
   KappaND(1,2) = std::complex<double>(pars(70), pars(71));
   KappaND(2,0) = std::complex<double>(pars(72), pars(73));
   KappaND(2,1) = std::complex<double>(pars(74), pars(75));
   KappaND(2,2) = std::complex<double>(pars(76), pars(77));

}

cSMHdCKM_soft_parameters::Soft_traces cSMHdCKM_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

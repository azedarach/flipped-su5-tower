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

// File generated at Mon 5 Nov 2018 12:48:08

#include "cSMHdCKMRHN_soft_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int cSMHdCKMRHN_soft_parameters::numberOfParameters;

cSMHdCKMRHN_soft_parameters::cSMHdCKMRHN_soft_parameters(const cSMHdCKMRHN_input_parameters& input_)
   : cSMHdCKMRHN_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

cSMHdCKMRHN_soft_parameters::cSMHdCKMRHN_soft_parameters(
   const cSMHdCKMRHN_susy_parameters& susy_model
   , const Eigen::Matrix<std::complex<double>,3,3>& Mv_, double mu2_, double v_
)
   : cSMHdCKMRHN_susy_parameters(susy_model)
   , Mv(Mv_), mu2(mu2_), v(v_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd cSMHdCKMRHN_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

cSMHdCKMRHN_soft_parameters cSMHdCKMRHN_soft_parameters::calc_beta(int loops) const
{
   Eigen::Matrix<std::complex<double>,3,3> beta_Mv = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   double beta_mu2 = 0.;
   double beta_v = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_Mv += calc_beta_Mv_1_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_Mv += calc_beta_Mv_2_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_2_loop(TRACE_STRUCT);

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


   const cSMHdCKMRHN_susy_parameters susy_betas(cSMHdCKMRHN_susy_parameters::calc_beta(loops));

   return cSMHdCKMRHN_soft_parameters(susy_betas, beta_Mv, beta_mu2, beta_v);
}

cSMHdCKMRHN_soft_parameters cSMHdCKMRHN_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void cSMHdCKMRHN_soft_parameters::clear()
{
   cSMHdCKMRHN_susy_parameters::clear();

   Mv = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   mu2 = 0.;
   v = 0.;

}

Eigen::ArrayXd cSMHdCKMRHN_soft_parameters::get() const
{
   Eigen::ArrayXd pars(cSMHdCKMRHN_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(76) = Re(Mv(0,0));
   pars(77) = Im(Mv(0,0));
   pars(78) = Re(Mv(0,1));
   pars(79) = Im(Mv(0,1));
   pars(80) = Re(Mv(0,2));
   pars(81) = Im(Mv(0,2));
   pars(82) = Re(Mv(1,0));
   pars(83) = Im(Mv(1,0));
   pars(84) = Re(Mv(1,1));
   pars(85) = Im(Mv(1,1));
   pars(86) = Re(Mv(1,2));
   pars(87) = Im(Mv(1,2));
   pars(88) = Re(Mv(2,0));
   pars(89) = Im(Mv(2,0));
   pars(90) = Re(Mv(2,1));
   pars(91) = Im(Mv(2,1));
   pars(92) = Re(Mv(2,2));
   pars(93) = Im(Mv(2,2));
   pars(94) = mu2;
   pars(95) = v;


   return pars;
}

void cSMHdCKMRHN_soft_parameters::print(std::ostream& ostr) const
{
   cSMHdCKMRHN_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "Mv = " << Mv << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';

}

void cSMHdCKMRHN_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   cSMHdCKMRHN_susy_parameters::set(pars);

   Mv(0,0) = std::complex<double>(pars(76), pars(77));
   Mv(0,1) = std::complex<double>(pars(78), pars(79));
   Mv(0,2) = std::complex<double>(pars(80), pars(81));
   Mv(1,0) = std::complex<double>(pars(82), pars(83));
   Mv(1,1) = std::complex<double>(pars(84), pars(85));
   Mv(1,2) = std::complex<double>(pars(86), pars(87));
   Mv(2,0) = std::complex<double>(pars(88), pars(89));
   Mv(2,1) = std::complex<double>(pars(90), pars(91));
   Mv(2,2) = std::complex<double>(pars(92), pars(93));
   mu2 = pars(94);
   v = pars(95);

}

cSMHdCKMRHN_soft_parameters::Soft_traces cSMHdCKMRHN_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYv = Re((Yv*Yv.adjoint()).trace());
      TRACE_STRUCT.traceMvconjMvYvAdjYv = Re((Mv*Mv.conjugate()*Yv*Yv.adjoint()).
         trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYvYvAdjYe = Re((Ye*Yv.adjoint()*Yv*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );
      TRACE_STRUCT.traceYvAdjYvYvAdjYv = Re((Yv*Yv.adjoint()*Yv*Yv.adjoint()).trace()
         );
      TRACE_STRUCT.traceMvconjMvYvAdjYeYeAdjYv = Re((Mv*Mv.conjugate()*Yv*Ye.adjoint(
         )*Ye*Yv.adjoint()).trace());
      TRACE_STRUCT.traceMvconjMvYvAdjYvYvAdjYv = Re((Mv*Mv.conjugate()*Yv*Yv.adjoint(
         )*Yv*Yv.adjoint()).trace());
      TRACE_STRUCT.traceMvconjYvTpYvconjMvYvAdjYv = Re((Mv*Yv.conjugate()*Yv.
         transpose()*Mv.conjugate()*Yv*Yv.adjoint()).trace());

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKMRHN_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

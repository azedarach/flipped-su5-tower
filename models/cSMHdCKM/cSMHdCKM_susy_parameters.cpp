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

// File generated at Mon 5 Nov 2018 12:47:12

#include "cSMHdCKM_susy_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME cSMHdCKM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces
#define TRACE_STRUCT_TYPE Susy_traces
#define CALCULATE_TRACES(l) calc_susy_traces(l);

const int cSMHdCKM_susy_parameters::numberOfParameters;

cSMHdCKM_susy_parameters::cSMHdCKM_susy_parameters(const cSMHdCKM_input_parameters& input_)
   : input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

cSMHdCKM_susy_parameters::cSMHdCKM_susy_parameters(
   double scale_, int loops_, int thresholds_,
   const cSMHdCKM_input_parameters& input_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<std
   ::complex<double>,3,3>& Yd_, const Eigen::Matrix<std::complex<double>,3,3>&
   Ye_, const Eigen::Matrix<std::complex<double>,3,3>& Yu_
)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambdax(Lambdax_), Yd(Yd_), Ye(Ye_), Yu(Yu_)
   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd cSMHdCKM_susy_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

cSMHdCKM_susy_parameters cSMHdCKM_susy_parameters::calc_beta(int loops) const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambdax = 0.;
   Eigen::Matrix<std::complex<double>,3,3> beta_Yd = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Ye = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();
   Eigen::Matrix<std::complex<double>,3,3> beta_Yu = Eigen::Matrix<std::complex<
      double>,3,3>::Zero();

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_g1 += calc_beta_g1_1_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_1_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_1_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_1_loop(TRACE_STRUCT);
      beta_Yd += calc_beta_Yd_1_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_1_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_g1 += calc_beta_g1_2_loop(TRACE_STRUCT);
         beta_g2 += calc_beta_g2_2_loop(TRACE_STRUCT);
         beta_g3 += calc_beta_g3_2_loop(TRACE_STRUCT);
         beta_Lambdax += calc_beta_Lambdax_2_loop(TRACE_STRUCT);
         beta_Yd += calc_beta_Yd_2_loop(TRACE_STRUCT);
         beta_Ye += calc_beta_Ye_2_loop(TRACE_STRUCT);
         beta_Yu += calc_beta_Yu_2_loop(TRACE_STRUCT);

         if (loops > 2) {
         #ifdef ENABLE_THREADS
            {
               auto fut_g1 = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_g1_3_loop(TRACE_STRUCT); });
               auto fut_g2 = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_g2_3_loop(TRACE_STRUCT); });
               auto fut_g3 = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_g3_3_loop(TRACE_STRUCT); });
               auto fut_Lambdax = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_Lambdax_3_loop(TRACE_STRUCT);
                  });
               auto fut_Yd = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_Yd_3_loop(TRACE_STRUCT); });
               auto fut_Ye = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_Ye_3_loop(TRACE_STRUCT); });
               auto fut_Yu = global_thread_pool().run_packaged_task([this, &
                  TRACE_STRUCT](){ return calc_beta_Yu_3_loop(TRACE_STRUCT); });

               beta_g1 += fut_g1.get();
               beta_g2 += fut_g2.get();
               beta_g3 += fut_g3.get();
               beta_Lambdax += fut_Lambdax.get();
               beta_Yd += fut_Yd.get();
               beta_Ye += fut_Ye.get();
               beta_Yu += fut_Yu.get();

            }
         #else
            beta_g1 += calc_beta_g1_3_loop(TRACE_STRUCT);
            beta_g2 += calc_beta_g2_3_loop(TRACE_STRUCT);
            beta_g3 += calc_beta_g3_3_loop(TRACE_STRUCT);
            beta_Lambdax += calc_beta_Lambdax_3_loop(TRACE_STRUCT);
            beta_Yd += calc_beta_Yd_3_loop(TRACE_STRUCT);
            beta_Ye += calc_beta_Ye_3_loop(TRACE_STRUCT);
            beta_Yu += calc_beta_Yu_3_loop(TRACE_STRUCT);
         #endif

            if (loops > 3) {
               beta_g3 += calc_beta_g3_4_loop(TRACE_STRUCT);
               beta_Lambdax += calc_beta_Lambdax_4_loop(TRACE_STRUCT);
               beta_Yu += calc_beta_Yu_4_loop(TRACE_STRUCT);

            }
         }
      }
   }


   return cSMHdCKM_susy_parameters(get_scale(), loops, get_thresholds(), input,
                    beta_g1, beta_g2, beta_g3, beta_Lambdax, beta_Yd, beta_Ye, beta_Yu);
}

cSMHdCKM_susy_parameters cSMHdCKM_susy_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void cSMHdCKM_susy_parameters::clear()
{
   reset();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambdax = 0.;
   Yd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ye = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Yu = Eigen::Matrix<std::complex<double>,3,3>::Zero();

}



Eigen::ArrayXd cSMHdCKM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambdax;
   pars(4) = Re(Yd(0,0));
   pars(5) = Im(Yd(0,0));
   pars(6) = Re(Yd(0,1));
   pars(7) = Im(Yd(0,1));
   pars(8) = Re(Yd(0,2));
   pars(9) = Im(Yd(0,2));
   pars(10) = Re(Yd(1,0));
   pars(11) = Im(Yd(1,0));
   pars(12) = Re(Yd(1,1));
   pars(13) = Im(Yd(1,1));
   pars(14) = Re(Yd(1,2));
   pars(15) = Im(Yd(1,2));
   pars(16) = Re(Yd(2,0));
   pars(17) = Im(Yd(2,0));
   pars(18) = Re(Yd(2,1));
   pars(19) = Im(Yd(2,1));
   pars(20) = Re(Yd(2,2));
   pars(21) = Im(Yd(2,2));
   pars(22) = Re(Ye(0,0));
   pars(23) = Im(Ye(0,0));
   pars(24) = Re(Ye(0,1));
   pars(25) = Im(Ye(0,1));
   pars(26) = Re(Ye(0,2));
   pars(27) = Im(Ye(0,2));
   pars(28) = Re(Ye(1,0));
   pars(29) = Im(Ye(1,0));
   pars(30) = Re(Ye(1,1));
   pars(31) = Im(Ye(1,1));
   pars(32) = Re(Ye(1,2));
   pars(33) = Im(Ye(1,2));
   pars(34) = Re(Ye(2,0));
   pars(35) = Im(Ye(2,0));
   pars(36) = Re(Ye(2,1));
   pars(37) = Im(Ye(2,1));
   pars(38) = Re(Ye(2,2));
   pars(39) = Im(Ye(2,2));
   pars(40) = Re(Yu(0,0));
   pars(41) = Im(Yu(0,0));
   pars(42) = Re(Yu(0,1));
   pars(43) = Im(Yu(0,1));
   pars(44) = Re(Yu(0,2));
   pars(45) = Im(Yu(0,2));
   pars(46) = Re(Yu(1,0));
   pars(47) = Im(Yu(1,0));
   pars(48) = Re(Yu(1,1));
   pars(49) = Im(Yu(1,1));
   pars(50) = Re(Yu(1,2));
   pars(51) = Im(Yu(1,2));
   pars(52) = Re(Yu(2,0));
   pars(53) = Im(Yu(2,0));
   pars(54) = Re(Yu(2,1));
   pars(55) = Im(Yu(2,1));
   pars(56) = Re(Yu(2,2));
   pars(57) = Im(Yu(2,2));


   return pars;
}

void cSMHdCKM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters at Q = " << get_scale() << ":\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Yu = " << Yu << '\n';

}

void cSMHdCKM_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambdax = pars(3);
   Yd(0,0) = std::complex<double>(pars(4), pars(5));
   Yd(0,1) = std::complex<double>(pars(6), pars(7));
   Yd(0,2) = std::complex<double>(pars(8), pars(9));
   Yd(1,0) = std::complex<double>(pars(10), pars(11));
   Yd(1,1) = std::complex<double>(pars(12), pars(13));
   Yd(1,2) = std::complex<double>(pars(14), pars(15));
   Yd(2,0) = std::complex<double>(pars(16), pars(17));
   Yd(2,1) = std::complex<double>(pars(18), pars(19));
   Yd(2,2) = std::complex<double>(pars(20), pars(21));
   Ye(0,0) = std::complex<double>(pars(22), pars(23));
   Ye(0,1) = std::complex<double>(pars(24), pars(25));
   Ye(0,2) = std::complex<double>(pars(26), pars(27));
   Ye(1,0) = std::complex<double>(pars(28), pars(29));
   Ye(1,1) = std::complex<double>(pars(30), pars(31));
   Ye(1,2) = std::complex<double>(pars(32), pars(33));
   Ye(2,0) = std::complex<double>(pars(34), pars(35));
   Ye(2,1) = std::complex<double>(pars(36), pars(37));
   Ye(2,2) = std::complex<double>(pars(38), pars(39));
   Yu(0,0) = std::complex<double>(pars(40), pars(41));
   Yu(0,1) = std::complex<double>(pars(42), pars(43));
   Yu(0,2) = std::complex<double>(pars(44), pars(45));
   Yu(1,0) = std::complex<double>(pars(46), pars(47));
   Yu(1,1) = std::complex<double>(pars(48), pars(49));
   Yu(1,2) = std::complex<double>(pars(50), pars(51));
   Yu(2,0) = std::complex<double>(pars(52), pars(53));
   Yu(2,1) = std::complex<double>(pars(54), pars(55));
   Yu(2,2) = std::complex<double>(pars(56), pars(57));

}

const cSMHdCKM_input_parameters& cSMHdCKM_susy_parameters::get_input() const
{
   return input;
}

cSMHdCKM_input_parameters& cSMHdCKM_susy_parameters::get_input()
{
   return input;
}

void cSMHdCKM_susy_parameters::set_input_parameters(const cSMHdCKM_input_parameters& input_)
{
   input = input_;
}

cSMHdCKM_susy_parameters::Susy_traces cSMHdCKM_susy_parameters::calc_susy_traces(int loops) const
{
   Susy_traces susy_traces;

   if (loops > 0) {
      

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()*
         Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd = Re((Yd*Yd.adjoint()*Yd*Yu.adjoint()*
         Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()*
         Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yu.adjoint()*
         Yu*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()*
         Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()*
         Yu*Yu.adjoint()).trace());

   }

   if (loops > 2) {

   }

   return susy_traces;
}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

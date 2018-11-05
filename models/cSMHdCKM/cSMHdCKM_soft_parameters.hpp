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

#ifndef cSMHdCKM_soft_parameters_H
#define cSMHdCKM_soft_parameters_H

#include "cSMHdCKM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class cSMHdCKM_soft_parameters : public cSMHdCKM_susy_parameters {
public:
   explicit cSMHdCKM_soft_parameters(const cSMHdCKM_input_parameters& input_ = cSMHdCKM_input_parameters());
   cSMHdCKM_soft_parameters(const cSMHdCKM_susy_parameters& , double mu2_, double v_,
                            const Eigen::Matrix<std::complex<double>,3,3>& Kappa_);
   cSMHdCKM_soft_parameters(const cSMHdCKM_soft_parameters&) = default;
   cSMHdCKM_soft_parameters(cSMHdCKM_soft_parameters&&) = default;
   virtual ~cSMHdCKM_soft_parameters() = default;
   cSMHdCKM_soft_parameters& operator=(const cSMHdCKM_soft_parameters&) = default;
   cSMHdCKM_soft_parameters& operator=(cSMHdCKM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   cSMHdCKM_soft_parameters calc_beta() const;
   cSMHdCKM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }
   void set_Kappa(const Eigen::Matrix<std::complex<double>,3,3>& Kappa_) { Kappa = Kappa_; }
   void set_Kappa(int i, int k, const std::complex<double>& value) { Kappa(i,k) = value; }

   double get_mu2() const { return mu2; }
   double get_v() const { return v; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Kappa() const { return Kappa; }
   std::complex<double> get_Kappa(int i, int k) const { return Kappa(i,k); }


protected:
   double mu2{};
   double v{};
   Eigen::Matrix<std::complex<double>,3,3> Kappa{};


private:
   static const int numberOfParameters = 78;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};

   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Kappa_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Kappa_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Kappa_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Kappa_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const cSMHdCKM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif

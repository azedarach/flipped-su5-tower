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

#ifndef cSMHdCKMRHN_soft_parameters_H
#define cSMHdCKMRHN_soft_parameters_H

#include "cSMHdCKMRHN_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class cSMHdCKMRHN_soft_parameters : public cSMHdCKMRHN_susy_parameters {
public:
   explicit cSMHdCKMRHN_soft_parameters(const cSMHdCKMRHN_input_parameters& input_ = cSMHdCKMRHN_input_parameters());
   cSMHdCKMRHN_soft_parameters(const cSMHdCKMRHN_susy_parameters& , const Eigen::Matrix<std::complex<double>,3,3>& Mv_, double mu2_, double v_);
   cSMHdCKMRHN_soft_parameters(const cSMHdCKMRHN_soft_parameters&) = default;
   cSMHdCKMRHN_soft_parameters(cSMHdCKMRHN_soft_parameters&&) = default;
   virtual ~cSMHdCKMRHN_soft_parameters() = default;
   cSMHdCKMRHN_soft_parameters& operator=(const cSMHdCKMRHN_soft_parameters&) = default;
   cSMHdCKMRHN_soft_parameters& operator=(cSMHdCKMRHN_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   cSMHdCKMRHN_soft_parameters calc_beta() const;
   cSMHdCKMRHN_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_Mv(const Eigen::Matrix<std::complex<double>,3,3>& Mv_) { Mv = Mv_; }
   void set_Mv(int i, int k, const std::complex<double>& value) { Mv(i,k) = value; }
   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   const Eigen::Matrix<std::complex<double>,3,3>& get_Mv() const { return Mv; }
   std::complex<double> get_Mv(int i, int k) const { return Mv(i,k); }
   double get_mu2() const { return mu2; }
   double get_v() const { return v; }


protected:
   Eigen::Matrix<std::complex<double>,3,3> Mv{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   double mu2{};
   double v{};


private:
   static const int numberOfParameters = 96;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYvAdjYv{};
      double traceMvconjMvYvAdjYv{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYeAdjYvYvAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double traceYvAdjYvYvAdjYv{};
      double traceMvconjMvYvAdjYeYeAdjYv{};
      double traceMvconjMvYvAdjYvYvAdjYv{};
      double traceMvconjYvTpYvconjMvYvAdjYv{};

   };
   Soft_traces calc_soft_traces(int) const;

   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Mv_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Mv_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Mv_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Mv_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const cSMHdCKMRHN_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif

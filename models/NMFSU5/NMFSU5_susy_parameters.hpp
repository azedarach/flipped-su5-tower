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

#ifndef NMFSU5_susy_parameters_H
#define NMFSU5_susy_parameters_H

#include "betafunction.hpp"
#include "NMFSU5_input_parameters.hpp"

#include <iosfwd>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Susy_traces

class NMFSU5_susy_parameters : public Beta_function {
public:
   explicit NMFSU5_susy_parameters(const NMFSU5_input_parameters& input_ = NMFSU5_input_parameters());
   NMFSU5_susy_parameters(double scale_, int loops_, int thresholds_, const NMFSU5_input_parameters& input_,
                          double g5_, double gX_, double lam1_, double lam2_, double lam3_, double lam3t_,
                          double lam4_, double lam4t_, double lam5_, double lam5t_, double lam6_, double lam6t_,
                          const std::complex<double>& lam7_, const std::complex<double>& lam8_,
                          const std::complex<double>& eta1_, const std::complex<double>& eta2_,
                          const std::complex<double>& eta3_, const Eigen::Matrix<std::complex<double>,3,3>& Y10_,
                          const Eigen::Matrix<std::complex<double>,3,3>& Y10Pr_,
                          const Eigen::Matrix<std::complex<double>,3,3>& Y5b_,
                          const Eigen::Matrix<std::complex<double>,3,3>& Y5bPr_,
                          const Eigen::Matrix<std::complex<double>,3,3>& Y1_,
                          const Eigen::Matrix<std::complex<double>,3,3>& Y1Pr_);
   NMFSU5_susy_parameters(const NMFSU5_susy_parameters&) = default;
   NMFSU5_susy_parameters(NMFSU5_susy_parameters&&) = default;
   virtual ~NMFSU5_susy_parameters() = default;
   NMFSU5_susy_parameters& operator=(const NMFSU5_susy_parameters&) = default;
   NMFSU5_susy_parameters& operator=(NMFSU5_susy_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&) override;
   const NMFSU5_input_parameters& get_input() const;
   NMFSU5_input_parameters& get_input();
   void set_input_parameters(const NMFSU5_input_parameters&);

   NMFSU5_susy_parameters calc_beta() const;
   NMFSU5_susy_parameters calc_beta(int) const;
   virtual void clear();

   void set_g5(double g5_) { g5 = g5_; }
   void set_gX(double gX_) { gX = gX_; }
   void set_lam1(double lam1_) { lam1 = lam1_; }
   void set_lam2(double lam2_) { lam2 = lam2_; }
   void set_lam3(double lam3_) { lam3 = lam3_; }
   void set_lam3t(double lam3t_) { lam3t = lam3t_; }
   void set_lam4(double lam4_) { lam4 = lam4_; }
   void set_lam4t(double lam4t_) { lam4t = lam4t_; }
   void set_lam5(double lam5_) { lam5 = lam5_; }
   void set_lam5t(double lam5t_) { lam5t = lam5t_; }
   void set_lam6(double lam6_) { lam6 = lam6_; }
   void set_lam6t(double lam6t_) { lam6t = lam6t_; }
   void set_lam7(const std::complex<double>& lam7_) { lam7 = lam7_; }
   void set_lam8(const std::complex<double>& lam8_) { lam8 = lam8_; }
   void set_eta1(const std::complex<double>& eta1_) { eta1 = eta1_; }
   void set_eta2(const std::complex<double>& eta2_) { eta2 = eta2_; }
   void set_eta3(const std::complex<double>& eta3_) { eta3 = eta3_; }
   void set_Y10(const Eigen::Matrix<std::complex<double>,3,3>& Y10_) { Y10 = Y10_; }
   void set_Y10(int i, int k, const std::complex<double>& value) { Y10(i,k) = value; }
   void set_Y10Pr(const Eigen::Matrix<std::complex<double>,3,3>& Y10Pr_) { Y10Pr = Y10Pr_; }
   void set_Y10Pr(int i, int k, const std::complex<double>& value) { Y10Pr(i,k) = value; }
   void set_Y5b(const Eigen::Matrix<std::complex<double>,3,3>& Y5b_) { Y5b = Y5b_; }
   void set_Y5b(int i, int k, const std::complex<double>& value) { Y5b(i,k) = value; }
   void set_Y5bPr(const Eigen::Matrix<std::complex<double>,3,3>& Y5bPr_) { Y5bPr = Y5bPr_; }
   void set_Y5bPr(int i, int k, const std::complex<double>& value) { Y5bPr(i,k) = value; }
   void set_Y1(const Eigen::Matrix<std::complex<double>,3,3>& Y1_) { Y1 = Y1_; }
   void set_Y1(int i, int k, const std::complex<double>& value) { Y1(i,k) = value; }
   void set_Y1Pr(const Eigen::Matrix<std::complex<double>,3,3>& Y1Pr_) { Y1Pr = Y1Pr_; }
   void set_Y1Pr(int i, int k, const std::complex<double>& value) { Y1Pr(i,k) = value; }

   double get_g5() const { return g5; }
   double get_gX() const { return gX; }
   double get_lam1() const { return lam1; }
   double get_lam2() const { return lam2; }
   double get_lam3() const { return lam3; }
   double get_lam3t() const { return lam3t; }
   double get_lam4() const { return lam4; }
   double get_lam4t() const { return lam4t; }
   double get_lam5() const { return lam5; }
   double get_lam5t() const { return lam5t; }
   double get_lam6() const { return lam6; }
   double get_lam6t() const { return lam6t; }
   std::complex<double> get_lam7() const { return lam7; }
   std::complex<double> get_lam8() const { return lam8; }
   std::complex<double> get_eta1() const { return eta1; }
   std::complex<double> get_eta2() const { return eta2; }
   std::complex<double> get_eta3() const { return eta3; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y10() const { return Y10; }
   std::complex<double> get_Y10(int i, int k) const { return Y10(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y10Pr() const { return Y10Pr; }
   std::complex<double> get_Y10Pr(int i, int k) const { return Y10Pr(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y5b() const { return Y5b; }
   std::complex<double> get_Y5b(int i, int k) const { return Y5b(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y5bPr() const { return Y5bPr; }
   std::complex<double> get_Y5bPr(int i, int k) const { return Y5bPr(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y1() const { return Y1; }
   std::complex<double> get_Y1(int i, int k) const { return Y1(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Y1Pr() const { return Y1Pr; }
   std::complex<double> get_Y1Pr(int i, int k) const { return Y1Pr(i,k); }


protected:
   double g5{};
   double gX{};
   double lam1{};
   double lam2{};
   double lam3{};
   double lam3t{};
   double lam4{};
   double lam4t{};
   double lam5{};
   double lam5t{};
   double lam6{};
   double lam6t{};
   std::complex<double> lam7{};
   std::complex<double> lam8{};
   std::complex<double> eta1{};
   std::complex<double> eta2{};
   std::complex<double> eta3{};
   Eigen::Matrix<std::complex<double>,3,3> Y10{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Y10Pr{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Y5b{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Y5bPr{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Y1{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Y1Pr{Eigen::Matrix<std::complex<double>,3,3>::Zero()};

   NMFSU5_input_parameters input{};

private:
   static const int numberOfParameters = 130;

   struct Susy_traces {
   };
   Susy_traces calc_susy_traces(int) const;

   double calc_beta_g5_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g5_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g5_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_g5_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gX_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gX_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gX_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_gX_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam1_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam1_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam1_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam1_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam2_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3t_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3t_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3t_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam3t_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4t_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4t_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4t_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam4t_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5t_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5t_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5t_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam5t_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6t_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6t_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6t_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_lam6t_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam7_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam7_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam7_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam7_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam8_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam8_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam8_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_lam8_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta1_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta1_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta1_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta1_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta2_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta2_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta2_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta2_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta3_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta3_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta3_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_eta3_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10Pr_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10Pr_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10Pr_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y10Pr_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5b_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5b_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5b_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5b_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5bPr_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5bPr_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5bPr_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y5bPr_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1_4_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1Pr_1_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1Pr_2_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1Pr_3_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<std::complex<double>,3,3> calc_beta_Y1Pr_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const NMFSU5_susy_parameters&);

#undef TRACE_STRUCT_TYPE

} // namespace flexiblesusy

#endif

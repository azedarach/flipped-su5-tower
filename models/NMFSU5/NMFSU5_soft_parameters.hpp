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

#ifndef NMFSU5_soft_parameters_H
#define NMFSU5_soft_parameters_H

#include "NMFSU5_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class NMFSU5_soft_parameters : public NMFSU5_susy_parameters {
public:
   explicit NMFSU5_soft_parameters(const NMFSU5_input_parameters& input_ = NMFSU5_input_parameters());
   NMFSU5_soft_parameters(const NMFSU5_susy_parameters& , double m10sq_, double m5sq_, double m5Prsq_,
                          const std::complex<double>& m12sq_, const std::complex<double>& mu_,
                          const std::complex<double>& muPr_, double v_, double vPr_, double VG_);
   NMFSU5_soft_parameters(const NMFSU5_soft_parameters&) = default;
   NMFSU5_soft_parameters(NMFSU5_soft_parameters&&) = default;
   virtual ~NMFSU5_soft_parameters() = default;
   NMFSU5_soft_parameters& operator=(const NMFSU5_soft_parameters&) = default;
   NMFSU5_soft_parameters& operator=(NMFSU5_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   NMFSU5_soft_parameters calc_beta() const;
   NMFSU5_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_m10sq(double m10sq_) { m10sq = m10sq_; }
   void set_m5sq(double m5sq_) { m5sq = m5sq_; }
   void set_m5Prsq(double m5Prsq_) { m5Prsq = m5Prsq_; }
   void set_m12sq(const std::complex<double>& m12sq_) { m12sq = m12sq_; }
   void set_mu(const std::complex<double>& mu_) { mu = mu_; }
   void set_muPr(const std::complex<double>& muPr_) { muPr = muPr_; }
   void set_v(double v_) { v = v_; }
   void set_vPr(double vPr_) { vPr = vPr_; }
   void set_VG(double VG_) { VG = VG_; }

   double get_m10sq() const { return m10sq; }
   double get_m5sq() const { return m5sq; }
   double get_m5Prsq() const { return m5Prsq; }
   std::complex<double> get_m12sq() const { return m12sq; }
   std::complex<double> get_mu() const { return mu; }
   std::complex<double> get_muPr() const { return muPr; }
   double get_v() const { return v; }
   double get_vPr() const { return vPr; }
   double get_VG() const { return VG; }

protected:
   double m10sq{};
   double m5sq{};
   double m5Prsq{};
   std::complex<double> m12sq{};
   std::complex<double> mu{};
   std::complex<double> muPr{};
   double v{};
   double vPr{};
   double VG{};

private:
   static const int numberOfParameters = 142;

   struct Soft_traces {
   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_m10sq_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m10sq_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m10sq_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m10sq_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5sq_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5sq_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5sq_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5sq_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5Prsq_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5Prsq_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5Prsq_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_m5Prsq_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_m12sq_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_m12sq_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_m12sq_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_m12sq_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_mu_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_mu_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_mu_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_mu_4_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_muPr_1_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_muPr_2_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_muPr_3_loop(const TRACE_STRUCT_TYPE&) const;
   std::complex<double> calc_beta_muPr_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vPr_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vPr_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vPr_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vPr_4_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_VG_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_VG_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_VG_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_VG_4_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const NMFSU5_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif

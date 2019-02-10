
#ifndef NMFSU5EFT_TWO_SCALE_LOW_SCALE_CONSTRAINT_H
#define NMFSU5EFT_TWO_SCALE_LOW_SCALE_CONSTRAINT_H

#include "NMFSU5EFT_low_scale_constraint.hpp"
#include "NMFSU5EFT_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class cSMHdCKM;

class Two_scale;

template<>
class NMFSU5EFT_low_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   NMFSU5EFT_low_scale_constraint() = default;
   NMFSU5EFT_low_scale_constraint(cSMHdCKM<Two_scale>*, const softsusy::QedQcd&,
                                  const NMFSU5EFT_input_parameters&);
   virtual ~NMFSU5EFT_low_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "NMFSU5 low-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm();
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns();
   double get_initial_scale_guess() const;
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);

private:
   double scale{0.};
   double initial_scale_guess{0.};
   cSMHdCKM<Two_scale>* model{nullptr};
   softsusy::QedQcd qedqcd{};
   NMFSU5EFT_input_parameters input{};
   Eigen::Matrix<std::complex<double>,3,3> ckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> pmns{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,3,3> neutrinoDRbar{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> neutrinoBasis{Eigen::Matrix<double,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> neutrinoMix{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   double mW_run{0.};
   double mZ_run{0.};
   double AlphaS{0.};
   double e_run{0.};
   double ThetaWDRbar{0.};
   double new_g1{0.}, new_g2{0.}, new_g3{0.};

   double calculate_theta_w();
   void calculate_threshold_corrections();
   void calculate_DRbar_gauge_couplings();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_Kappa();
   void calculate_running_SM_masses();
   void calculate_neutrino_basis();
   void calculate_neutrino_mixing();
   double calculate_delta_alpha_em(double) const;
   double calculate_delta_alpha_s(double) const;
   double calculate_alpha_s_SM5_at(softsusy::QedQcd, double) const;
   void check_model_ptr() const;
   void update_scale();
};

} // namespace flexiblesusy

#endif

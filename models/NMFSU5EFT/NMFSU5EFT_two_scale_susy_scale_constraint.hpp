#ifndef NMFSU5EFT_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define NMFSU5EFT_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "NMFSU5EFT_susy_scale_constraint.hpp"
#include "NMFSU5EFT_input_parameters.hpp"
#include "single_scale_constraint.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <Eigen/LU>

#include <complex>

namespace flexiblesusy {

template <class T>
class cSMHdCKMRHN;

class Two_scale;

template<>
class NMFSU5EFT_susy_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   NMFSU5EFT_susy_scale_constraint() = default;
   NMFSU5EFT_susy_scale_constraint(cSMHdCKMRHN<Two_scale>*, const softsusy::QedQcd&,
      const NMFSU5EFT_input_parameters&);
   virtual ~NMFSU5EFT_susy_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override { return "NMFSU5EFT SUSY-scale constraint"; }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const;
   const NMFSU5EFT_input_parameters& get_input_parameters() const;
   cSMHdCKMRHN<Two_scale>* get_model() const;
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);

protected:
   void update_scale();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   cSMHdCKMRHN<Two_scale>* model{nullptr};
   softsusy::QedQcd qedqcd{};
   NMFSU5EFT_input_parameters input{};

   double calculate_initial_scale_guess() const;
   void calculate_seesaw_Mv(const Eigen::Matrix<std::complex<double>,3,3>&,
                            const Eigen::Matrix<std::complex<double>,3,3>&,
                            Eigen::Matrix<std::complex<double>,3,3>&) const;
   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif

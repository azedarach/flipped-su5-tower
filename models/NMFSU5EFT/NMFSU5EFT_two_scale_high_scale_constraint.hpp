#ifndef NMFSU5EFT_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define NMFSU5EFT_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "NMFSU5EFT_high_scale_constraint.hpp"
#include "NMFSU5EFT_input_parameters.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class cSMHdCKMRHN;

class Two_scale;

template<>
class NMFSU5EFT_high_scale_constraint<Two_scale> : public Single_scale_constraint {
public:
   NMFSU5EFT_high_scale_constraint() = default;
   NMFSU5EFT_high_scale_constraint(cSMHdCKMRHN<Two_scale>*,
                                   const NMFSU5EFT_input_parameters&);
   virtual ~NMFSU5EFT_high_scale_constraint() = default;
   virtual void apply() override;
   virtual double get_scale() const override;
   virtual std::string name() const override {
      return "NMFSU5EFT high-scale constraint";
   }
   virtual void set_model(Model*) override;

   void clear();
   double get_initial_scale_guess() const;
   const NMFSU5EFT_input_parameters& get_input_parameters() const;
   cSMHdCKMRHN<Two_scale>* get_model() const;
   void initialize();
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();
   bool check_non_perturbative();
   bool check_high_scale_non_perturbative();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   cSMHdCKMRHN<Two_scale>* model{nullptr};
   NMFSU5<Two_scale>* high_scale_model{nullptr};
   NMFSU5EFT_input_parameters input{};

   void check_model_ptr() const;
};

} // namespace flexiblesusy

#endif

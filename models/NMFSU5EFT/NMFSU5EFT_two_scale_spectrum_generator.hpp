#ifndef NMFSU5EFT_TWO_SCALE_SPECTRUM_GENERATOR_H
#define NMFSU5EFT_TWO_SCALE_SPECTRUM_GENERATOR_H

#include "NMFSU5EFT_spectrum_generator.hpp"
#include "NMFSU5EFT_spectrum_generator_interface.hpp"

#include "NMFSU5_two_scale_model.hpp"
#include "NMFSU5_model_slha.hpp"

#include "cSMHdCKMRHN_two_scale_model.hpp"
#include "cSMHdCKMRHN_model_slha.hpp"

#include "cSMHdCKM_two_scale_model.hpp"
#include "cSMHdCKM_model_slha.hpp"

namespace softsusy { class QedQcd; }

namespace flexiblesusy {

class Two_scale;

template <>
class NMFSU5EFT_spectrum_generator<Two_scale>
   : public NMFSU5EFT_spectrum_generator_interface<Two_scale> {
public:
   NMFSU5EFT_spectrum_generator() = default;
   virtual ~NMFSU5EFT_spectrum_generator() = default;

   double get_high_scale() const { return high_scale; }
   double get_susy_scale() const { return susy_scale; }
   double get_low_scale() const { return low_scale; }
   double get_pole_mass_scale() const { return get_pole_mass_scale(susy_scale); }

   void write_running_couplings(
      const std::string& filename = "NMFSU5EFT_rgflow.dat") const;

protected:
   virtual void run_except(
      const softsusy::QedQcd&,
      const NMFSU5EFT_input_parameters&) override;

private:
   double high_scale{0.};
   double susy_scale{0.};
   double low_scale{0.};

   void calculate_spectrum(double, double);
   double get_eft_pole_mass_scale(double, double) const;
   double get_pole_mass_scale(double) const;
};

} // namespace flexiblesusy

#endif

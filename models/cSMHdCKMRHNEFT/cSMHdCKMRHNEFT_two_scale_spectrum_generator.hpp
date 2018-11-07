#ifndef cSMHdCKMRHNEFT_TWO_SCALE_SPECTRUM_GENERATOR_H
#define cSMHdCKMRHNEFT_TWO_SCALE_SPECTRUM_GENERATOR_H

#include "cSMHdCKMRHNEFT_spectrum_generator.hpp"
#include "cSMHdCKMRHNEFT_spectrum_generator_interface.hpp"

#include "cSMHdCKMRHN_two_scale_model.hpp"

#include "cSMHdCKM_two_scale_model.hpp"

namespace softsusy { class QedQcd; }

namespace flexiblesusy {

class Two_scale;

template <>
class cSMHdCKMRHNEFT_spectrum_generator<Two_scale>
   : public cSMHdCKMRHNEFT_spectrum_generator_interface<Two_scale> {
public:
   cSMHdCKMRHNEFT_spectrum_generator() = default;
   virtual ~cSMHdCKMRHNEFT_spectrum_generator() = default;

   double get_high_scale() const { return high_scale; }
   double get_susy_scale() const { return susy_scale; }
   double get_low_scale() const { return low_scale; }
   double get_pole_mass_scale() const { return get_pole_mass_scale(susy_scale); }

   void write_running_couplings(const std::string& filename = "cSMHdCKMRHNEFT_rgflow.dat") const;

protected:
   virtual void run_except(
      const softsusy::QedQcd&, const cSMHdCKM_input_parameters&,
      const cSMHdCKMRHN_input_parameters&) override;

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

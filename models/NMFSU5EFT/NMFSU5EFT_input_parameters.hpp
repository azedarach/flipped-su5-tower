#ifndef NMFSU5EFT_INPUT_PARAMETERS_H
#define NMFSU5EFT_INPUT_PARAMETERS_H

#include "wrappers.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NMFSU5EFT_input_parameters {
   double Lambda1IN{};
   double Lambda2IN{};
   double Lambda3IN{};
   double Lambda3tIN{};
   double Lambda4IN{};
   double Lambda4tIN{};
   double Lambda5IN{};
   double Lambda5tIN{};
   double Lambda6IN{};
   double Lambda6tIN{};
   std::complex<double> Lambda7IN{};
   std::complex<double> Lambda8IN{};
   std::complex<double> Eta1IN{};
   std::complex<double> Eta2IN{};
   std::complex<double> Eta3IN{};
   int sign_delta_mAsq{1};

   Eigen::Matrix<std::complex<double>, 3, 3> UvInput{
      Eigen::Matrix<std::complex<double>, 3, 3>::Zero()};

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const NMFSU5EFT_input_parameters&);

} // namespace flexiblesusy

#endif

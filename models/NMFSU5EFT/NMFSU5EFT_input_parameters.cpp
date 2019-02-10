#include "NMFSU5EFT_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd NMFSU5EFT_input_parameters::get() const
{
   Eigen::ArrayXd pars(39);

   pars(0) = Lambda1IN;
   pars(1) = Lambda2IN;
   pars(2) = Lambda3IN;
   pars(3) = Lambda3tIN;
   pars(4) = Lambda4IN;
   pars(5) = Lambda4tIN;
   pars(6) = Lambda5IN;
   pars(7) = Lambda5tIN;
   pars(8) = Lambda6IN;
   pars(9) = Lambda6tIN;
   pars(10) = Re(Lambda7IN);
   pars(11) = Im(Lambda7IN);
   pars(12) = Re(Lambda8IN);
   pars(13) = Im(Lambda8IN);
   pars(14) = Re(Eta1IN);
   pars(15) = Im(Eta1IN);
   pars(16) = Re(Eta2IN);
   pars(17) = Im(Eta2IN);
   pars(18) = Re(Eta3IN);
   pars(19) = Im(Eta3IN);
   pars(20) = sign_delta_mAsq;
   pars(21) = Re(UvInput(0,0));
   pars(22) = Im(UvInput(0,0));
   pars(23) = Re(UvInput(0,1));
   pars(24) = Im(UvInput(0,1));
   pars(25) = Re(UvInput(0,2));
   pars(26) = Im(UvInput(0,2));
   pars(27) = Re(UvInput(1,0));
   pars(28) = Im(UvInput(1,0));
   pars(29) = Re(UvInput(1,1));
   pars(30) = Im(UvInput(1,1));
   pars(31) = Re(UvInput(1,2));
   pars(32) = Im(UvInput(1,2));
   pars(33) = Re(UvInput(2,0));
   pars(34) = Im(UvInput(2,0));
   pars(35) = Re(UvInput(2,1));
   pars(36) = Im(UvInput(2,1));
   pars(37) = Re(UvInput(2,2));
   pars(38) = Im(UvInput(2,2));

   return pars;
}

void NMFSU5EFT_input_parameters::set(const Eigen::ArrayXd& pars)
{
   Lambda1IN = pars(0);
   Lambda2IN = pars(1);
   Lambda3IN = pars(2);
   Lambda3tIN = pars(3);
   Lambda4IN = pars(4);
   Lambda4tIN = pars(5);
   Lambda5IN = pars(6);
   Lambda5tIN = pars(7);
   Lambda6IN = pars(8);
   Lambda6tIN = pars(8);
   Lambda7IN = std::complex<double>(pars(10), pars(11));
   Lambda8IN = std::complex<double>(pars(12), pars(13));
   Eta1IN = std::complex<double>(pars(14), pars(15));
   Eta2IN = std::complex<double>(pars(16), pars(17));
   Eta3IN = std::complex<double>(pars(18), pars(19));
   sign_delta_mAsq = pars(20);
   UvInput(0,0) = std::complex<double>(pars(21), pars(22));
   UvInput(0,1) = std::complex<double>(pars(23), pars(24));
   UvInput(0,2) = std::complex<double>(pars(25), pars(26));
   UvInput(1,0) = std::complex<double>(pars(27), pars(28));
   UvInput(1,1) = std::complex<double>(pars(29), pars(30));
   UvInput(1,2) = std::complex<double>(pars(31), pars(32));
   UvInput(2,0) = std::complex<double>(pars(33), pars(34));
   UvInput(2,1) = std::complex<double>(pars(35), pars(36));
   UvInput(2,2) = std::complex<double>(pars(37), pars(38));
}

std::ostream& operator<<(std::ostream& ostr, const NMFSU5EFT_input_parameters& input)
{
   ostr << "Lambda1IN = " << INPUT(Lambda1IN) << ", ";
   ostr << "Lambda2IN = " << INPUT(Lambda2IN) << ", ";
   ostr << "Lambda3IN = " << INPUT(Lambda3IN) << ", ";
   ostr << "Lambda3tIN = " << INPUT(Lambda3tIN) << ", ";
   ostr << "Lambda4IN = " << INPUT(Lambda4IN) << ", ";
   ostr << "Lambda4tIN = " << INPUT(Lambda4tIN) << ", ";
   ostr << "Lambda5IN = " << INPUT(Lambda5IN) << ", ";
   ostr << "Lambda5tIN = " << INPUT(Lambda5tIN) << ", ";
   ostr << "Lambda6IN = " << INPUT(Lambda6IN) << ", ";
   ostr << "Lambda6tIN = " << INPUT(Lambda6tIN) << ", ";
   ostr << "Lambda7IN = " << INPUT(Lambda7IN) << ", ";
   ostr << "Lambda8IN = " << INPUT(Lambda8IN) << ", ";
   ostr << "Eta1IN = " << INPUT(Eta1IN) << ", ";
   ostr << "Eta2IN = " << INPUT(Eta2IN) << ", ";
   ostr << "Eta3IN = " << INPUT(Eta3IN) << ", ";
   ostr << "sign_delta_mAsq = " << INPUT(sign_delta_mAsq) << ", ";
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         ostr << "UvInput(" << i << ", " << j << ") = " << INPUT(UvInput(i, j)) << ", ";
      }
   }

   return ostr;
}

} // namespace flexiblesusy

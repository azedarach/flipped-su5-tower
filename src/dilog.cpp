// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "dilog.hpp"
#include "numerics2.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {
   template <typename T>
   T sqr(T x) noexcept { return x*x; }
} // namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 */
double dilog(double x) noexcept {
   const double PI = M_PI;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 *
 * @note not full long double precision yet!
 */
long double dilog(long double x) noexcept {
   return dilog(static_cast<double>(x));
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<long double> dilog(const std::complex<long double>& z) noexcept {
   using flexiblesusy::fast_log;
   const long double PI = 3.1415926535897932384626433832795l;
   std::complex<long double> cy, cz;
   int jsgn, ipi12;
   static const int N = 20;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 19}]
   const long double bf[N] = {
      - 1.l/4.l,
      + 1.l/36.l,
      - 1.l/36.e2l,
      + 1.l/21168.e1l,
      - 1.l/108864.e2l,
      + 1.l/52690176.e1l,
      - 4.0647616451442255268059093862919666745470571274397078e-11l,
      + 8.9216910204564525552179873167527488515142836130490451e-13l,
      - 1.9939295860721075687236443477937897056306947496538801e-14l,
      + 4.5189800296199181916504765528555932283968190144666184e-16l,
      - 1.0356517612181247014483411542218656665960912381686505e-17l,
      + 2.3952186210261867457402837430009803816789490019429743e-19l,
      - 5.5817858743250093362830745056254199055670546676443981e-21l,
      + 1.3091507554183212858123073991865923017498498387833038e-22l,
      - 3.0874198024267402932422797648664624315955652561327457e-24l,
      + 7.315975652702203420357905609252148591033401063690875e-26l,
      - 1.7408456572340007409890551477597025453408414217542713e-27l,
      + 4.1576356446138997196178996207752266734882541595115639e-29l,
      - 9.9621484882846221031940067024558388498548600173944888e-31l,
      + 2.3940344248961653005211679878937495629342791569329158e-32l,
   };

   const long double rz = std::real(z);
   const long double iz = std::imag(z);
   const long double az = std::sqrt(sqr(rz) + sqr(iz));

   // special cases
   if (iz == 0.l) {
      if (rz <= 1.l)
         return std::complex<long double>(dilog(rz), 0.l);
      if (rz > 1.l)
         return std::complex<long double>(dilog(rz), -PI*std::log(rz));
   } else if (az < std::numeric_limits<long double>::epsilon()) {
      return z;
   }

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5l) {
      if (az > 1.l) {
         cy = -0.5l * sqr(fast_log(-z));
         cz = -fast_log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -fast_log(1.l - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= std::sqrt(2*rz)) {
         cz = -fast_log(z);
         cy = cz * fast_log(1.l - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5l * sqr(fast_log(-z));
         cz = -fast_log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<long double> cz2(sqr(cz));
   std::complex<long double> sumC;

   for (int i1 = 2; i1 < N; i1++)
      sumC = cz2 * (sumC + bf[N + 1 - i1]);

   // lowest order terms w/ different powers
   sumC = cz + cz2 * (bf[0] + cz * (bf[1] + sumC));

   const std::complex<long double> result
      = static_cast<long double>(jsgn) * sumC + cy + ipi12 * PI * PI / 12.l;

   return result;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<double> dilog(const std::complex<double>& z) noexcept {
   const auto rz = static_cast<long double>(std::real(z));
   const auto iz = static_cast<long double>(std::imag(z));
   const auto re = dilog(std::complex<long double>(rz, iz));
   return std::complex<double>(std::real(re), std::imag(re));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta)\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
double clausen_2(double x) noexcept
{
   using std::exp;
   const std::complex<double> img(0.,1.);

   return std::imag(dilog(exp(img*x)));
}

} // namespace flexiblesusy

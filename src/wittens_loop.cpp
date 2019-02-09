#include "wittens_loop.hpp"
#include "dilog.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

namespace {

constexpr double piSqrOverSix = Pi * Pi / 6.0;

double fbb(double b) noexcept
{
   const auto rad = Sqrt(1. - 4. * b);
   const auto lb = Log(b);
   const auto lrad = Log(-(rad - 1.) / (rad + 1.));

   return -((2. * b - 1.) * (2. * dilog((rad - 1.) / (rad + 1.)) + piSqrOverSix
                             + 0.5 * lrad * lrad) / rad) - 0.5 * lb * lb;
}

double fzerob(double b) noexcept
{
   return dilog(1. - b);
}

double fzerobinv(double b) noexcept
{
   const auto lb = Log(b);
   return -0.5 * lb * lb - fzerob(b);
}

double f(double a, double b) noexcept
{
   if (is_zero(a) && is_zero(b)) {
      return piSqrOverSix;
   } else if (is_zero(a)) {
      return fzerob(b);
   } else if (is_equal(a, b)) {
      return fbb(a);
   }

   const auto q = 1. - 2. * (a + b) + (a - b) * (a - b);
   const auto x1 = 0.5 * (1. + b - a + Sqrt(q));
   const auto x2 = 0.5 * (1. + b - a - Sqrt(q));
   const auto y1 = 0.5 * (1. + a - b + Sqrt(q));
   const auto y2 = 0.5 * (1. + a - b - Sqrt(q));

   const auto sum = dilog(-x2 / y1) + dilog(-y2 / x1) - dilog(-x1 / y2)
      - dilog(-y1 / x2) + dilog((b - a) / x2) + dilog((a - b)/ y2)
      - dilog((b - a) / x1) - dilog((a - b) / y1);

   return -0.5 * Log(a) * Log(b) + 0.5 * (1. - a - b) * sum / Sqrt(q);
}

} // anonymous namespace

double calc_I3(double mDelta_sq, double mX_sq) noexcept
{
   const auto s = mDelta_sq / mX_sq;
   const auto s_inv = mX_sq / mDelta_sq;
   const auto ls = Log(s);
   const auto f00 = piSqrOverSix;

   return 1. + 2. * ls + s * (1. - 2. * s) * ls * ls
      + 2. * (s_inv - 1.) * (
         f00 * (1. + s + s * s) + 2. * s * f(1., s)
         + fzerob(s) * (1. + s) * (1. + 2. * s) + s * s * f(s_inv, s_inv));
}

} // namespace flexiblesusy

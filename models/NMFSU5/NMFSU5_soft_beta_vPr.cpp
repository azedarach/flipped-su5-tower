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

#include "NMFSU5_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of vPr.
 *
 * @return 1-loop beta function
 */
double NMFSU5_soft_parameters::calc_beta_vPr_1_loop(const Soft_traces& /* soft_traces */) const
{


   double beta_vPr;

   beta_vPr = 0;


   return beta_vPr;
}

/**
 * Calculates the 2-loop beta function of vPr.
 *
 * @return 2-loop beta function
 */
double NMFSU5_soft_parameters::calc_beta_vPr_2_loop(const Soft_traces& /* soft_traces */) const
{


   double beta_vPr;

   beta_vPr = 0;


   return beta_vPr;
}

/**
 * Calculates the 3-loop beta function of vPr.
 *
 * @return 3-loop beta function
 */
double NMFSU5_soft_parameters::calc_beta_vPr_3_loop(const Soft_traces& /* soft_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vPr;

   beta_vPr = 0;


   return beta_vPr;
}

/**
 * Calculates the 4-loop beta function of vPr.
 *
 * @return 4-loop beta function
 */
double NMFSU5_soft_parameters::calc_beta_vPr_4_loop(const Soft_traces& /* soft_traces */) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vPr;

   beta_vPr = 0;


   return beta_vPr;
}

} // namespace flexiblesusy
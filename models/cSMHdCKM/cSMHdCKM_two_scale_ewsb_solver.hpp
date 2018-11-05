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

// File generated at Mon 5 Nov 2018 12:47:54

/**
 * @file cSMHdCKM_two_scale_ewsb_solver.hpp
 *
 * @brief contains class for solving EWSB when two-scale algorithm is used
 *
 * This file was generated at Mon 5 Nov 2018 12:47:54 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#ifndef cSMHdCKM_TWO_SCALE_EWSB_SOLVER_H
#define cSMHdCKM_TWO_SCALE_EWSB_SOLVER_H

#include "cSMHdCKM_ewsb_solver.hpp"
#include "cSMHdCKM_ewsb_solver_interface.hpp"
#include "error.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Two_scale;

class cSMHdCKM_mass_eigenstates;

template<>
class cSMHdCKM_ewsb_solver<Two_scale> : public cSMHdCKM_ewsb_solver_interface {
public:
   cSMHdCKM_ewsb_solver() = default;
   cSMHdCKM_ewsb_solver(const cSMHdCKM_ewsb_solver&) = default;
   cSMHdCKM_ewsb_solver(cSMHdCKM_ewsb_solver&&) = default;
   virtual ~cSMHdCKM_ewsb_solver() {}
   cSMHdCKM_ewsb_solver& operator=(const cSMHdCKM_ewsb_solver&) = default;
   cSMHdCKM_ewsb_solver& operator=(cSMHdCKM_ewsb_solver&&) = default;

   virtual void set_loop_order(int l) override { loop_order = l; }
   virtual void set_number_of_iterations(int n) override { number_of_iterations = n; }
   virtual void set_precision(double p) override { precision = p; }

   virtual int get_loop_order() const override { return loop_order; }
   virtual int get_number_of_iterations() const override { return number_of_iterations; }
   virtual double get_precision() const override { return precision; }

   virtual int solve(cSMHdCKM_mass_eigenstates&) override;
private:
   static const int number_of_ewsb_equations = 1;
   using EWSB_vector_t = Eigen::Matrix<double,number_of_ewsb_equations,1>;

   class EEWSBStepFailed : public Error {
   public:
      virtual ~EEWSBStepFailed() {}
      virtual std::string what() const { return "Could not perform EWSB step."; }
   };

   int number_of_iterations{100}; ///< maximum number of iterations
   int loop_order{2};             ///< loop order to solve EWSB at
   double precision{1.e-5};       ///< precision goal

   void set_ewsb_solution(cSMHdCKM_mass_eigenstates&, const EWSB_solver*);
   template <typename It> void set_best_ewsb_solution(cSMHdCKM_mass_eigenstates&, It, It);

   int solve_tree_level(cSMHdCKM_mass_eigenstates&);
   int solve_iteratively(cSMHdCKM_mass_eigenstates&);
   int solve_iteratively_at(cSMHdCKM_mass_eigenstates&, int);
   int solve_iteratively_with(cSMHdCKM_mass_eigenstates&, EWSB_solver*, const EWSB_vector_t&);

   EWSB_vector_t initial_guess(const cSMHdCKM_mass_eigenstates&) const;
   EWSB_vector_t tadpole_equations(const cSMHdCKM_mass_eigenstates&) const;
   EWSB_vector_t ewsb_step(const cSMHdCKM_mass_eigenstates&) const;
};

} // namespace flexiblesusy

#endif

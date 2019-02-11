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

/**
 * @file NMFSU5_mass_eigenstates.cpp
 * @brief implementation of the NMFSU5 model class
 *
 * Contains the definition of the NMFSU5 model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated at Mon 5 Nov 2018 12:47:51 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#include "NMFSU5_mass_eigenstates.hpp"
#include "NMFSU5_ewsb_solver_interface.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "ewsb_solver.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "pv.hpp"
#include "raii.hpp"
#include "thread_pool.hpp"
#include "functors.hpp"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "NMFSU5_two_scale_ewsb_solver.hpp"
#endif




#include "sm_threeloop_as.hpp"


#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>
#include <stdexcept>

namespace flexiblesusy {

#define CLASSNAME NMFSU5_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS       loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS       loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT       loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU   loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION            loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS    loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS    loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_SCHEME                 loop_corrections.higgs_3L_scheme
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS    loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT    loop_corrections.higgs_at_at_at
#define HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS loop_corrections.higgs_at_as_as_as

CLASSNAME::NMFSU5_mass_eigenstates(const NMFSU5_input_parameters& input_)
   : NMFSU5_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new NMFSU5_ewsb_solver<Two_scale>())
#endif
{
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::do_calculate_bsm_pole_masses(bool flag)
{
   calculate_bsm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_bsm_pole_masses() const
{
   return calculate_bsm_pole_masses;
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

void CLASSNAME::set_ewsb_loop_order(int loop_order)
{
   ewsb_loop_order = loop_order;
   if (ewsb_solver) {
      ewsb_solver->set_loop_order(ewsb_loop_order);
   }
}

void CLASSNAME::set_loop_corrections(const Loop_corrections& loop_corrections_)
{
   loop_corrections = loop_corrections_;
}

const Loop_corrections& CLASSNAME::get_loop_corrections() const
{
   return loop_corrections;
}

void CLASSNAME::set_threshold_corrections(const Threshold_corrections& tc)
{
   threshold_corrections = tc;
}

const Threshold_corrections& CLASSNAME::get_threshold_corrections() const
{
   return threshold_corrections;
}

int CLASSNAME::get_number_of_ewsb_iterations() const
{
   return static_cast<int>(std::abs(-log10(ewsb_iteration_precision) * 10));
}

int CLASSNAME::get_number_of_mass_iterations() const
{
   return static_cast<int>(std::abs(-log10(precision) * 10));
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision_);
   }
}

void CLASSNAME::set_pole_mass_loop_order(int loop_order)
{
   pole_mass_loop_order = loop_order;
}

int CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision);
   }
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const NMFSU5_physical& CLASSNAME::get_physical() const
{
   return physical;
}

NMFSU5_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const NMFSU5_physical& physical_)
{
   physical = physical_;
}

const Problems& CLASSNAME::get_problems() const
{
   return problems;
}

Problems& CLASSNAME::get_problems()
{
   return problems;
}

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<NMFSU5_ewsb_solver_interface>& solver)
{
   ewsb_solver = solver;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   const auto tadpole_(tadpole_equations());
   std::copy(tadpole_.data(), tadpole_.data() + number_of_ewsb_equations, tadpole);
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole(
      Eigen::Matrix<double,number_of_ewsb_equations,1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop(0));
      tadpole[1] -= Re(tadpole_hh_1loop(1));
      tadpole[2] -= Re(tadpole_hh_1loop(2));

      if (ewsb_loop_order > 1) {

      }
   }

   return tadpole;
}

/**
 * This function returns the vector of tadpoles, each divided by the
 * corresponding VEV.  Thus, the returned tadpoles have the dimension
 * GeV^2 each.
 *
 * @return vector of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations_over_vevs() const
{
   auto tadpole = tadpole_equations();

   tadpole[0] /= VG;
   tadpole[1] /= v;
   tadpole[2] /= vPr;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;



   return error;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   if (!ewsb_solver) {
      throw SetupError("NMFSU5_mass_eigenstates::solve_ewsb_tree_level: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(0);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb_one_loop()
{
   if (!ewsb_solver) {
      throw SetupError("NMFSU5_mass_eigenstates::solve_ewsb_one_loop: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(1);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb()
{
   if (!ewsb_solver) {
      throw SetupError("NMFSU5_mass_eigenstates::solve_ewsb: "
                       "no EWSB solver set");
   }

   VERBOSE_MSG("\t\tSolving NMFSU5 EWSB at " << ewsb_loop_order << "-loop order");

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(ewsb_loop_order);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "NMFSU5\n"
           "========================================\n";
   NMFSU5_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "UV = " << UV << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 * @note: They take squared arguments!
 */

double CLASSNAME::A0(double m) const noexcept
{
   return passarino_veltman::ReA0(m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB1(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB00(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB22(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReH0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReF0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReG0(p, m1, m2, Sqr(get_scale()));
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
//   const auto save_mu2_raii = make_raii_save(mu2);

   const bool has_no_ewsb_flag = problems.no_ewsb();
   const auto save_ewsb_flag = make_raii_guard(
      [this, has_no_ewsb_flag] () {
         if (has_no_ewsb_flag) {
            this->problems.flag_no_ewsb_tree_level();
         } else {
            this->problems.unflag_no_ewsb_tree_level();
         }
      }
   );
   problems.unflag_no_ewsb_tree_level();
   solve_ewsb_tree_level();
#ifdef ENABLE_VERBOSE
   if (problems.no_ewsb()) {
      WARNING("solving EWSB at 0-loop order failed");
   }
#endif

   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MFv();

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 9u));

   if (calculate_bsm_pole_masses) {
   }

   if (calculate_sm_pole_masses) {
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
   }

   if (calculate_sm_pole_masses) {
      calculate_MFv_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(UV) = UV;

}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{

}

/**
 * calculates spectrum for model once the MSbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MFv = Eigen::Matrix<double,6,1>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   UV = Eigen::Matrix<std::complex<double>,6,6>::Zero();


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   NMFSU5_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MFv(0) = pars(0);
   MFv(1) = pars(1);
   MFv(2) = pars(2);
   MFv(3) = pars(3);
   MFv(4) = pars(4);
   MFv(5) = pars(5);
   MFd(0) = pars(6);
   MFd(1) = pars(7);
   MFd(2) = pars(8);
   MFu(0) = pars(9);
   MFu(1) = pars(10);
   MFu(2) = pars(11);
   MFe(0) = pars(12);
   MFe(1) = pars(13);
   MFe(2) = pars(14);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(15);

   pars(0) = MFv(0);
   pars(1) = MFv(1);
   pars(2) = MFv(2);
   pars(3) = MFv(3);
   pars(4) = MFv(4);
   pars(5) = MFv(5);
   pars(6) = MFd(0);
   pars(7) = MFd(1);
   pars(8) = MFd(2);
   pars(9) = MFu(0);
   pars(10) = MFu(1);
   pars(11) = MFu(2);
   pars(12) = MFe(0);
   pars(13) = MFe(1);
   pars(14) = MFe(2);

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   Vd(0,0) = std::complex<double>(pars(15), pars(16));
   Vd(0,1) = std::complex<double>(pars(17), pars(18));
   Vd(0,2) = std::complex<double>(pars(19), pars(20));
   Vd(1,0) = std::complex<double>(pars(21), pars(22));
   Vd(1,1) = std::complex<double>(pars(23), pars(24));
   Vd(1,2) = std::complex<double>(pars(25), pars(26));
   Vd(2,0) = std::complex<double>(pars(27), pars(28));
   Vd(2,1) = std::complex<double>(pars(29), pars(30));
   Vd(2,2) = std::complex<double>(pars(31), pars(32));
   Ud(0,0) = std::complex<double>(pars(33), pars(34));
   Ud(0,1) = std::complex<double>(pars(35), pars(36));
   Ud(0,2) = std::complex<double>(pars(37), pars(38));
   Ud(1,0) = std::complex<double>(pars(39), pars(40));
   Ud(1,1) = std::complex<double>(pars(41), pars(42));
   Ud(1,2) = std::complex<double>(pars(43), pars(44));
   Ud(2,0) = std::complex<double>(pars(45), pars(46));
   Ud(2,1) = std::complex<double>(pars(47), pars(48));
   Ud(2,2) = std::complex<double>(pars(49), pars(50));
   Vu(0,0) = std::complex<double>(pars(51), pars(52));
   Vu(0,1) = std::complex<double>(pars(53), pars(54));
   Vu(0,2) = std::complex<double>(pars(55), pars(56));
   Vu(1,0) = std::complex<double>(pars(57), pars(58));
   Vu(1,1) = std::complex<double>(pars(59), pars(60));
   Vu(1,2) = std::complex<double>(pars(61), pars(62));
   Vu(2,0) = std::complex<double>(pars(63), pars(64));
   Vu(2,1) = std::complex<double>(pars(65), pars(66));
   Vu(2,2) = std::complex<double>(pars(67), pars(68));
   Uu(0,0) = std::complex<double>(pars(69), pars(70));
   Uu(0,1) = std::complex<double>(pars(71), pars(72));
   Uu(0,2) = std::complex<double>(pars(73), pars(74));
   Uu(1,0) = std::complex<double>(pars(75), pars(76));
   Uu(1,1) = std::complex<double>(pars(77), pars(78));
   Uu(1,2) = std::complex<double>(pars(79), pars(80));
   Uu(2,0) = std::complex<double>(pars(81), pars(82));
   Uu(2,1) = std::complex<double>(pars(83), pars(84));
   Uu(2,2) = std::complex<double>(pars(85), pars(86));
   Ve(0,0) = std::complex<double>(pars(87), pars(88));
   Ve(0,1) = std::complex<double>(pars(89), pars(90));
   Ve(0,2) = std::complex<double>(pars(91), pars(92));
   Ve(1,0) = std::complex<double>(pars(93), pars(94));
   Ve(1,1) = std::complex<double>(pars(95), pars(96));
   Ve(1,2) = std::complex<double>(pars(97), pars(98));
   Ve(2,0) = std::complex<double>(pars(99), pars(100));
   Ve(2,1) = std::complex<double>(pars(101), pars(102));
   Ve(2,2) = std::complex<double>(pars(103), pars(104));
   Ue(0,0) = std::complex<double>(pars(105), pars(106));
   Ue(0,1) = std::complex<double>(pars(107), pars(108));
   Ue(0,2) = std::complex<double>(pars(109), pars(110));
   Ue(1,0) = std::complex<double>(pars(111), pars(112));
   Ue(1,1) = std::complex<double>(pars(113), pars(114));
   Ue(1,2) = std::complex<double>(pars(115), pars(116));
   Ue(2,0) = std::complex<double>(pars(117), pars(118));
   Ue(2,1) = std::complex<double>(pars(119), pars(120));
   Ue(2,2) = std::complex<double>(pars(121), pars(122));
   UV(0,0) = std::complex<double>(pars(123), pars(124));
   UV(0,1) = std::complex<double>(pars(125), pars(126));
   UV(0,2) = std::complex<double>(pars(127), pars(128));
   UV(0,3) = std::complex<double>(pars(129), pars(130));
   UV(0,4) = std::complex<double>(pars(131), pars(132));
   UV(0,5) = std::complex<double>(pars(133), pars(134));
   UV(1,0) = std::complex<double>(pars(135), pars(136));
   UV(1,1) = std::complex<double>(pars(137), pars(138));
   UV(1,2) = std::complex<double>(pars(139), pars(140));
   UV(1,3) = std::complex<double>(pars(141), pars(142));
   UV(1,4) = std::complex<double>(pars(143), pars(144));
   UV(1,5) = std::complex<double>(pars(145), pars(146));
   UV(2,0) = std::complex<double>(pars(147), pars(148));
   UV(2,1) = std::complex<double>(pars(149), pars(150));
   UV(2,2) = std::complex<double>(pars(151), pars(152));
   UV(2,3) = std::complex<double>(pars(153), pars(154));
   UV(2,4) = std::complex<double>(pars(155), pars(156));
   UV(2,5) = std::complex<double>(pars(157), pars(158));
   UV(3,0) = std::complex<double>(pars(159), pars(160));
   UV(3,1) = std::complex<double>(pars(161), pars(162));
   UV(3,2) = std::complex<double>(pars(163), pars(164));
   UV(3,3) = std::complex<double>(pars(165), pars(166));
   UV(3,4) = std::complex<double>(pars(167), pars(168));
   UV(3,5) = std::complex<double>(pars(169), pars(170));
   UV(4,0) = std::complex<double>(pars(171), pars(172));
   UV(4,1) = std::complex<double>(pars(173), pars(174));
   UV(4,2) = std::complex<double>(pars(175), pars(176));
   UV(4,3) = std::complex<double>(pars(177), pars(178));
   UV(4,4) = std::complex<double>(pars(179), pars(180));
   UV(4,5) = std::complex<double>(pars(181), pars(182));
   UV(5,0) = std::complex<double>(pars(183), pars(184));
   UV(5,1) = std::complex<double>(pars(185), pars(186));
   UV(5,2) = std::complex<double>(pars(187), pars(188));
   UV(5,3) = std::complex<double>(pars(189), pars(190));
   UV(5,4) = std::complex<double>(pars(191), pars(192));
   UV(5,5) = std::complex<double>(pars(193), pars(194));

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(195);

   pars(15) = Re(Vd(0,0));
   pars(16) = Im(Vd(0,0));
   pars(17) = Re(Vd(0,1));
   pars(18) = Im(Vd(0,1));
   pars(19) = Re(Vd(0,2));
   pars(20) = Im(Vd(0,2));
   pars(21) = Re(Vd(1,0));
   pars(22) = Im(Vd(1,0));
   pars(23) = Re(Vd(1,1));
   pars(24) = Im(Vd(1,1));
   pars(25) = Re(Vd(1,2));
   pars(26) = Im(Vd(1,2));
   pars(27) = Re(Vd(2,0));
   pars(28) = Im(Vd(2,0));
   pars(29) = Re(Vd(2,1));
   pars(30) = Im(Vd(2,1));
   pars(31) = Re(Vd(2,2));
   pars(32) = Im(Vd(2,2));
   pars(33) = Re(Ud(0,0));
   pars(34) = Im(Ud(0,0));
   pars(35) = Re(Ud(0,1));
   pars(36) = Im(Ud(0,1));
   pars(37) = Re(Ud(0,2));
   pars(38) = Im(Ud(0,2));
   pars(39) = Re(Ud(1,0));
   pars(40) = Im(Ud(1,0));
   pars(41) = Re(Ud(1,1));
   pars(42) = Im(Ud(1,1));
   pars(43) = Re(Ud(1,2));
   pars(44) = Im(Ud(1,2));
   pars(45) = Re(Ud(2,0));
   pars(46) = Im(Ud(2,0));
   pars(47) = Re(Ud(2,1));
   pars(48) = Im(Ud(2,1));
   pars(49) = Re(Ud(2,2));
   pars(50) = Im(Ud(2,2));
   pars(51) = Re(Vu(0,0));
   pars(52) = Im(Vu(0,0));
   pars(53) = Re(Vu(0,1));
   pars(54) = Im(Vu(0,1));
   pars(55) = Re(Vu(0,2));
   pars(56) = Im(Vu(0,2));
   pars(57) = Re(Vu(1,0));
   pars(58) = Im(Vu(1,0));
   pars(59) = Re(Vu(1,1));
   pars(60) = Im(Vu(1,1));
   pars(61) = Re(Vu(1,2));
   pars(62) = Im(Vu(1,2));
   pars(63) = Re(Vu(2,0));
   pars(64) = Im(Vu(2,0));
   pars(65) = Re(Vu(2,1));
   pars(66) = Im(Vu(2,1));
   pars(67) = Re(Vu(2,2));
   pars(68) = Im(Vu(2,2));
   pars(69) = Re(Uu(0,0));
   pars(70) = Im(Uu(0,0));
   pars(71) = Re(Uu(0,1));
   pars(72) = Im(Uu(0,1));
   pars(73) = Re(Uu(0,2));
   pars(74) = Im(Uu(0,2));
   pars(75) = Re(Uu(1,0));
   pars(76) = Im(Uu(1,0));
   pars(77) = Re(Uu(1,1));
   pars(78) = Im(Uu(1,1));
   pars(79) = Re(Uu(1,2));
   pars(80) = Im(Uu(1,2));
   pars(81) = Re(Uu(2,0));
   pars(82) = Im(Uu(2,0));
   pars(83) = Re(Uu(2,1));
   pars(84) = Im(Uu(2,1));
   pars(85) = Re(Uu(2,2));
   pars(86) = Im(Uu(2,2));
   pars(87) = Re(Ve(0,0));
   pars(88) = Im(Ve(0,0));
   pars(89) = Re(Ve(0,1));
   pars(90) = Im(Ve(0,1));
   pars(91) = Re(Ve(0,2));
   pars(92) = Im(Ve(0,2));
   pars(93) = Re(Ve(1,0));
   pars(94) = Im(Ve(1,0));
   pars(95) = Re(Ve(1,1));
   pars(96) = Im(Ve(1,1));
   pars(97) = Re(Ve(1,2));
   pars(98) = Im(Ve(1,2));
   pars(99) = Re(Ve(2,0));
   pars(100) = Im(Ve(2,0));
   pars(101) = Re(Ve(2,1));
   pars(102) = Im(Ve(2,1));
   pars(103) = Re(Ve(2,2));
   pars(104) = Im(Ve(2,2));
   pars(105) = Re(Ue(0,0));
   pars(106) = Im(Ue(0,0));
   pars(107) = Re(Ue(0,1));
   pars(108) = Im(Ue(0,1));
   pars(109) = Re(Ue(0,2));
   pars(110) = Im(Ue(0,2));
   pars(111) = Re(Ue(1,0));
   pars(112) = Im(Ue(1,0));
   pars(113) = Re(Ue(1,1));
   pars(114) = Im(Ue(1,1));
   pars(115) = Re(Ue(1,2));
   pars(116) = Im(Ue(1,2));
   pars(117) = Re(Ue(2,0));
   pars(118) = Im(Ue(2,0));
   pars(119) = Re(Ue(2,1));
   pars(120) = Im(Ue(2,1));
   pars(121) = Re(Ue(2,2));
   pars(122) = Im(Ue(2,2));
   pars(123) = Re(UV(0,0));
   pars(124) = Im(UV(0,0));
   pars(125) = Re(UV(0,1));
   pars(126) = Im(UV(0,1));
   pars(127) = Re(UV(0,2));
   pars(128) = Im(UV(0,2));
   pars(129) = Re(UV(0,3));
   pars(130) = Im(UV(0,3));
   pars(131) = Re(UV(0,4));
   pars(132) = Im(UV(0,4));
   pars(133) = Re(UV(0,5));
   pars(134) = Im(UV(0,5));
   pars(135) = Re(UV(1,0));
   pars(136) = Im(UV(1,0));
   pars(137) = Re(UV(1,1));
   pars(138) = Im(UV(1,1));
   pars(139) = Re(UV(1,2));
   pars(140) = Im(UV(1,2));
   pars(141) = Re(UV(1,3));
   pars(142) = Im(UV(1,3));
   pars(143) = Re(UV(1,4));
   pars(144) = Im(UV(1,4));
   pars(145) = Re(UV(1,5));
   pars(146) = Im(UV(1,5));
   pars(147) = Re(UV(2,0));
   pars(148) = Im(UV(2,0));
   pars(149) = Re(UV(2,1));
   pars(150) = Im(UV(2,1));
   pars(151) = Re(UV(2,2));
   pars(152) = Im(UV(2,2));
   pars(153) = Re(UV(2,3));
   pars(154) = Im(UV(2,3));
   pars(155) = Re(UV(2,4));
   pars(156) = Im(UV(2,4));
   pars(157) = Re(UV(2,5));
   pars(158) = Im(UV(2,5));
   pars(159) = Re(UV(3,0));
   pars(160) = Im(UV(3,0));
   pars(161) = Re(UV(3,1));
   pars(162) = Im(UV(3,1));
   pars(163) = Re(UV(3,2));
   pars(164) = Im(UV(3,2));
   pars(165) = Re(UV(3,3));
   pars(166) = Im(UV(3,3));
   pars(167) = Re(UV(3,4));
   pars(168) = Im(UV(3,4));
   pars(169) = Re(UV(3,5));
   pars(170) = Im(UV(3,5));
   pars(171) = Re(UV(4,0));
   pars(172) = Im(UV(4,0));
   pars(173) = Re(UV(4,1));
   pars(174) = Im(UV(4,1));
   pars(175) = Re(UV(4,2));
   pars(176) = Im(UV(4,2));
   pars(177) = Re(UV(4,3));
   pars(178) = Im(UV(4,3));
   pars(179) = Re(UV(4,4));
   pars(180) = Im(UV(4,4));
   pars(181) = Re(UV(4,5));
   pars(182) = Im(UV(4,5));
   pars(183) = Re(UV(5,0));
   pars(184) = Im(UV(5,0));
   pars(185) = Re(UV(5,1));
   pars(186) = Im(UV(5,1));
   pars(187) = Re(UV(5,2));
   pars(188) = Im(UV(5,2));
   pars(189) = Re(UV(5,3));
   pars(190) = Im(UV(5,3));
   pars(191) = Re(UV(5,4));
   pars(192) = Im(UV(5,4));
   pars(193) = Re(UV(5,5));
   pars(194) = Im(UV(5,5));

   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

}

std::string CLASSNAME::name() const
{
   return "NMFSU5";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   NMFSU5_soft_parameters::run_to(scale, eps);
}








Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<std::complex<double>,6,6> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(0,3) = -0.7071067811865475 * (Y5b(0,0) * v + Y5bPr(0,0) * vPr);
   mass_matrix_Fv(0,4) = -0.7071067811865475 * (Y5b(1,0) * v + Y5bPr(1,0) * vPr);
   mass_matrix_Fv(0,5) = -0.7071067811865475 * (Y5b(2,0) * v + Y5bPr(2,0) * vPr);
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(1,3) = -0.7071067811865475 * (Y5b(0,1) * v + Y5bPr(0,1) * vPr);
   mass_matrix_Fv(1,4) = -0.7071067811865475 * (Y5b(1,1) * v + Y5bPr(1,1) * vPr);
   mass_matrix_Fv(1,5) = -0.7071067811865475 * (Y5b(2,1) * v + Y5bPr(2,1) * vPr);
   mass_matrix_Fv(2,2) = 0;
   mass_matrix_Fv(2,3) = -0.7071067811865475 * (Y5b(0,2) * v + Y5bPr(0,2) * vPr);
   mass_matrix_Fv(2,4) = -0.7071067811865475 * (Y5b(1,2) * v + Y5bPr(1,1) * vPr);
   mass_matrix_Fv(2,5) = -0.7071067811865475 * (Y5b(2,2) * v + Y5bPr(2,2) * vPr);
   mass_matrix_Fv(3,3) = 0;
   mass_matrix_Fv(3,4) = 0;
   mass_matrix_Fv(3,5) = 0;
   mass_matrix_Fv(4,4) = 0;
   mass_matrix_Fv(4,5) = 0;
   mass_matrix_Fv(5,5) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   const auto mass_matrix_Fv(get_mass_matrix_Fv());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV, eigenvalue_error);
   problems.flag_bad_mass(NMFSU5_info::Fv, eigenvalue_error > precision * Abs(
                             MFv(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV);
#endif
   normalize_to_interval(UV);

}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = -0.5 * (Y10(0,0) * v + Y10Pr(0,0) * vPr);
   mass_matrix_Fd(0,1) = -0.25 * (Y10(0,1) * v + Y10(1,0) * v
                                  + Y10Pr(0,1) * vPr + Y10Pr(1,0) * vPr);
   mass_matrix_Fd(0,2) = -0.25 * (Y10(0,2) * v + Y10(2,0) * v
                                  + Y10Pr(0,2) * vPr + Y10Pr(2,0) * vPr);
   mass_matrix_Fd(1,1) = -0.5 * (Y10(1,1) * v + Y10Pr(1,1) * vPr);
   mass_matrix_Fd(1,2) = -0.25 * (Y10(1,2) * v + Y10(2,1) * v
                                  + Y10Pr(1,2) * vPr + Y10Pr(2,1) * vPr);
   mass_matrix_Fd(2,2) = -0.5 * (Y10(2,2) * v + Y10Pr(2,2) * vPr);

   Symmetrize(mass_matrix_Fd);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(NMFSU5_info::Fd, eigenvalue_error > precision * Abs
      (MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475 * (Y5b(0,0) * v + Y5bPr(0,0) * vPr);
   mass_matrix_Fu(0,1) = -0.7071067811865475 * (Y5b(0,1) * v + Y5bPr(0,1) * vPr);
   mass_matrix_Fu(0,2) = -0.7071067811865475 * (Y5b(0,2) * v + Y5bPr(0,2) * vPr);
   mass_matrix_Fu(1,0) = -0.7071067811865475 * (Y5b(1,0) * v + Y5bPr(1,0) * vPr);
   mass_matrix_Fu(1,1) = -0.7071067811865475 * (Y5b(1,1) * v + Y5bPr(1,1) * vPr);
   mass_matrix_Fu(1,2) = -0.7071067811865475 * (Y5b(1,2) * v + Y5bPr(1,2) * vPr);
   mass_matrix_Fu(2,0) = -0.7071067811865475 * (Y5b(2,0) * v + Y5bPr(2,0) * vPr);
   mass_matrix_Fu(2,1) = -0.7071067811865475 * (Y5b(2,1) * v + Y5bPr(2,1) * vPr);
   mass_matrix_Fu(2,2) = -0.7071067811865475 * (Y5b(2,2) * v + Y5bPr(2,2) * vPr);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(NMFSU5_info::Fu, eigenvalue_error > precision * Abs
      (MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = Y1(0,0) * v + Y1Pr(0,0) * vPr;
   mass_matrix_Fe(0,1) = Y1(0,1) * v + Y1Pr(0,1) * vPr;
   mass_matrix_Fe(0,2) = Y1(0,2) * v + Y1Pr(0,2) * vPr;
   mass_matrix_Fe(1,0) = Y1(1,0) * v + Y1Pr(1,0) * vPr;
   mass_matrix_Fe(1,1) = Y1(1,1) * v + Y1Pr(1,1) * vPr;
   mass_matrix_Fe(1,2) = Y1(1,2) * v + Y1Pr(1,2) * vPr;
   mass_matrix_Fe(2,0) = Y1(2,0) * v + Y1Pr(2,0) * vPr;
   mass_matrix_Fe(2,1) = Y1(2,1) * v + Y1Pr(2,1) * vPr;
   mass_matrix_Fe(2,2) = Y1(2,2) * v + Y1Pr(2,2) * vPr;

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(NMFSU5_info::Fe, eigenvalue_error > precision * Abs
      (MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = 0;

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = 0;

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = 0;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Fv_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Fv_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Fv_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result = 0;

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Fv_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PL(p, i, k);

   return self_energy;
}


std::complex<double> CLASSNAME::tadpole_hh_1loop(int gO1) const
{
   std::complex<double> result = 0;


   return result * oneOver16PiSqr;
}




void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<std::complex<double>,6,6> M_tree(get_mass_matrix_Fv());
   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MFv(es));
      const Eigen::Matrix<std::complex<double>,6,6> self_energy_1  =
         self_energy_Fv_1loop_1(p);
      const Eigen::Matrix<std::complex<double>,6,6> self_energy_PL =
         self_energy_Fv_1loop_PL(p);
      const Eigen::Matrix<std::complex<double>,6,6> self_energy_PR =
         self_energy_Fv_1loop_PR(p);
      const Eigen::Matrix<std::complex<double>,6,6> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<std::complex<double>,6,6> M_loop(M_tree + 0.5 * (
         delta_M + delta_M.transpose()));
      Eigen::Array<double,6,1> eigen_values;
      decltype(UV) mix_UV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_UV,
            eigenvalue_error);
         problems.flag_bad_mass(cSMHdCKMRHN_info::Fv, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else

         fs_diagonalize_symmetric(M_loop, eigen_values, mix_UV);
      #endif
         normalize_to_interval(mix_UV);
      if (es == 0)
         PHYSICAL(UV) = mix_UV;
      PHYSICAL(MFv(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<std::complex<double>,3,3> M_tree(get_mass_matrix_Fd());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_1  =
         self_energy_Fd_1loop_1(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PL =
         self_energy_Fd_1loop_PL(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PR =
         self_energy_Fd_1loop_PR(p);
      const Eigen::Matrix<std::complex<double>,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<std::complex<double>,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(NMFSU5_info::Fd, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<std::complex<double>,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_1  =
         self_energy_Fu_1loop_1(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PL =
         self_energy_Fu_1loop_PL(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PR =
         self_energy_Fu_1loop_PR(p);
      const Eigen::Matrix<std::complex<double>,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<std::complex<double>,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Vu;
      decltype(Ue) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(NMFSU5_info::Fu, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Vu;
         PHYSICAL(Ue) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<std::complex<double>,3,3> M_tree(get_mass_matrix_Fe());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_1  =
         self_energy_Fe_1loop_1(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PL =
         self_energy_Fe_1loop_PL(p);
      const Eigen::Matrix<std::complex<double>,3,3> self_energy_PR =
         self_energy_Fe_1loop_PR(p);
      const Eigen::Matrix<std::complex<double>,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<std::complex<double>,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(NMFSU5_info::Fe, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}



double CLASSNAME::calculate_MFv_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1 = Re(self_energy_Fv_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fv_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fv_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL
                                                             + self_energy_PR);

   return m_drbar;
}

double CLASSNAME::calculate_MFe_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL
                                                             + self_energy_PR);

   return m_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL
                                                             + self_energy_PR);

   return m_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL
                                                             + self_energy_PR);

   return m_drbar;
}


std::ostream& operator<<(std::ostream& ostr, const NMFSU5_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

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

// File generated at Mon 5 Nov 2018 12:48:52

/**
 * @file cSMHdCKMRHN_mass_eigenstates.cpp
 * @brief implementation of the cSMHdCKMRHN model class
 *
 * Contains the definition of the cSMHdCKMRHN model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated at Mon 5 Nov 2018 12:48:52 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#include "cSMHdCKMRHN_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_ewsb_solver_interface.hpp"
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
#include "cSMHdCKMRHN_two_scale_ewsb_solver.hpp"
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

#define CLASSNAME cSMHdCKMRHN_mass_eigenstates

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

CLASSNAME::cSMHdCKMRHN_mass_eigenstates(const cSMHdCKMRHN_input_parameters& input_)
   : cSMHdCKMRHN_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new cSMHdCKMRHN_ewsb_solver<Two_scale>())
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

const cSMHdCKMRHN_physical& CLASSNAME::get_physical() const
{
   return physical;
}

cSMHdCKMRHN_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const cSMHdCKMRHN_physical& physical_)
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

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<cSMHdCKMRHN_ewsb_solver_interface>& solver)
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

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop());

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

   tadpole[0] /= v;


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
      throw SetupError("cSMHdCKMRHN_mass_eigenstates::solve_ewsb_tree_level: "
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
      throw SetupError("cSMHdCKMRHN_mass_eigenstates::solve_ewsb_one_loop: "
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
      throw SetupError("cSMHdCKMRHN_mass_eigenstates::solve_ewsb: "
                       "no EWSB solver set");
   }

   VERBOSE_MSG("\t\tSolving cSMHdCKMRHN EWSB at " << ewsb_loop_order << "-loop order");

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
           "cSMHdCKMRHN\n"
           "========================================\n";
   cSMHdCKMRHN_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHm = " << MHm << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

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
   ostr << "ZZ = " << ZZ << '\n';

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
   const auto save_mu2_raii = make_raii_save(mu2);

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

   calculate_MVPVZ();
   calculate_MVWp();
   calculate_MFv();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHm();
   calculate_MVG();

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
      tp.run_task([this] () { calculate_MVG_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MVP_pole(); });
      tp.run_task([this] () { calculate_MVZ_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_MVWp_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_Mhh_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MFv_pole();
      calculate_MVWp_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHm) = MHm;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(UV) = UV;
   PHYSICAL(MVWp) = MVWp;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;

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
   if (PHYSICAL(Mhh) < 0.) problems.flag_pole_tachyon(cSMHdCKMRHN_info::hh);

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
   MVG = 0.;
   MHm = 0.;
   MAh = 0.;
   Mhh = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFv = Eigen::Matrix<double,6,1>::Zero();
   UV = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MVWp = 0.;
   MVP = 0.;
   MVZ = 0.;



}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   cSMHdCKMRHN_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MHm = pars(1);
   MAh = pars(2);
   Mhh = pars(3);
   MFd(0) = pars(4);
   MFd(1) = pars(5);
   MFd(2) = pars(6);
   MFu(0) = pars(7);
   MFu(1) = pars(8);
   MFu(2) = pars(9);
   MFe(0) = pars(10);
   MFe(1) = pars(11);
   MFe(2) = pars(12);
   MFv(0) = pars(13);
   MFv(1) = pars(14);
   MFv(2) = pars(15);
   MFv(3) = pars(16);
   MFv(4) = pars(17);
   MFv(5) = pars(18);
   MVWp = pars(19);
   MVP = pars(20);
   MVZ = pars(21);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(22);

   pars(0) = MVG;
   pars(1) = MHm;
   pars(2) = MAh;
   pars(3) = Mhh;
   pars(4) = MFd(0);
   pars(5) = MFd(1);
   pars(6) = MFd(2);
   pars(7) = MFu(0);
   pars(8) = MFu(1);
   pars(9) = MFu(2);
   pars(10) = MFe(0);
   pars(11) = MFe(1);
   pars(12) = MFe(2);
   pars(13) = MFv(0);
   pars(14) = MFv(1);
   pars(15) = MFv(2);
   pars(16) = MFv(3);
   pars(17) = MFv(4);
   pars(18) = MFv(5);
   pars(19) = MVWp;
   pars(20) = MVP;
   pars(21) = MVZ;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   Vd(0,0) = std::complex<double>(pars(22), pars(23));
   Vd(0,1) = std::complex<double>(pars(24), pars(25));
   Vd(0,2) = std::complex<double>(pars(26), pars(27));
   Vd(1,0) = std::complex<double>(pars(28), pars(29));
   Vd(1,1) = std::complex<double>(pars(30), pars(31));
   Vd(1,2) = std::complex<double>(pars(32), pars(33));
   Vd(2,0) = std::complex<double>(pars(34), pars(35));
   Vd(2,1) = std::complex<double>(pars(36), pars(37));
   Vd(2,2) = std::complex<double>(pars(38), pars(39));
   Ud(0,0) = std::complex<double>(pars(40), pars(41));
   Ud(0,1) = std::complex<double>(pars(42), pars(43));
   Ud(0,2) = std::complex<double>(pars(44), pars(45));
   Ud(1,0) = std::complex<double>(pars(46), pars(47));
   Ud(1,1) = std::complex<double>(pars(48), pars(49));
   Ud(1,2) = std::complex<double>(pars(50), pars(51));
   Ud(2,0) = std::complex<double>(pars(52), pars(53));
   Ud(2,1) = std::complex<double>(pars(54), pars(55));
   Ud(2,2) = std::complex<double>(pars(56), pars(57));
   Vu(0,0) = std::complex<double>(pars(58), pars(59));
   Vu(0,1) = std::complex<double>(pars(60), pars(61));
   Vu(0,2) = std::complex<double>(pars(62), pars(63));
   Vu(1,0) = std::complex<double>(pars(64), pars(65));
   Vu(1,1) = std::complex<double>(pars(66), pars(67));
   Vu(1,2) = std::complex<double>(pars(68), pars(69));
   Vu(2,0) = std::complex<double>(pars(70), pars(71));
   Vu(2,1) = std::complex<double>(pars(72), pars(73));
   Vu(2,2) = std::complex<double>(pars(74), pars(75));
   Uu(0,0) = std::complex<double>(pars(76), pars(77));
   Uu(0,1) = std::complex<double>(pars(78), pars(79));
   Uu(0,2) = std::complex<double>(pars(80), pars(81));
   Uu(1,0) = std::complex<double>(pars(82), pars(83));
   Uu(1,1) = std::complex<double>(pars(84), pars(85));
   Uu(1,2) = std::complex<double>(pars(86), pars(87));
   Uu(2,0) = std::complex<double>(pars(88), pars(89));
   Uu(2,1) = std::complex<double>(pars(90), pars(91));
   Uu(2,2) = std::complex<double>(pars(92), pars(93));
   Ve(0,0) = std::complex<double>(pars(94), pars(95));
   Ve(0,1) = std::complex<double>(pars(96), pars(97));
   Ve(0,2) = std::complex<double>(pars(98), pars(99));
   Ve(1,0) = std::complex<double>(pars(100), pars(101));
   Ve(1,1) = std::complex<double>(pars(102), pars(103));
   Ve(1,2) = std::complex<double>(pars(104), pars(105));
   Ve(2,0) = std::complex<double>(pars(106), pars(107));
   Ve(2,1) = std::complex<double>(pars(108), pars(109));
   Ve(2,2) = std::complex<double>(pars(110), pars(111));
   Ue(0,0) = std::complex<double>(pars(112), pars(113));
   Ue(0,1) = std::complex<double>(pars(114), pars(115));
   Ue(0,2) = std::complex<double>(pars(116), pars(117));
   Ue(1,0) = std::complex<double>(pars(118), pars(119));
   Ue(1,1) = std::complex<double>(pars(120), pars(121));
   Ue(1,2) = std::complex<double>(pars(122), pars(123));
   Ue(2,0) = std::complex<double>(pars(124), pars(125));
   Ue(2,1) = std::complex<double>(pars(126), pars(127));
   Ue(2,2) = std::complex<double>(pars(128), pars(129));
   UV(0,0) = std::complex<double>(pars(130), pars(131));
   UV(0,1) = std::complex<double>(pars(132), pars(133));
   UV(0,2) = std::complex<double>(pars(134), pars(135));
   UV(0,3) = std::complex<double>(pars(136), pars(137));
   UV(0,4) = std::complex<double>(pars(138), pars(139));
   UV(0,5) = std::complex<double>(pars(140), pars(141));
   UV(1,0) = std::complex<double>(pars(142), pars(143));
   UV(1,1) = std::complex<double>(pars(144), pars(145));
   UV(1,2) = std::complex<double>(pars(146), pars(147));
   UV(1,3) = std::complex<double>(pars(148), pars(149));
   UV(1,4) = std::complex<double>(pars(150), pars(151));
   UV(1,5) = std::complex<double>(pars(152), pars(153));
   UV(2,0) = std::complex<double>(pars(154), pars(155));
   UV(2,1) = std::complex<double>(pars(156), pars(157));
   UV(2,2) = std::complex<double>(pars(158), pars(159));
   UV(2,3) = std::complex<double>(pars(160), pars(161));
   UV(2,4) = std::complex<double>(pars(162), pars(163));
   UV(2,5) = std::complex<double>(pars(164), pars(165));
   UV(3,0) = std::complex<double>(pars(166), pars(167));
   UV(3,1) = std::complex<double>(pars(168), pars(169));
   UV(3,2) = std::complex<double>(pars(170), pars(171));
   UV(3,3) = std::complex<double>(pars(172), pars(173));
   UV(3,4) = std::complex<double>(pars(174), pars(175));
   UV(3,5) = std::complex<double>(pars(176), pars(177));
   UV(4,0) = std::complex<double>(pars(178), pars(179));
   UV(4,1) = std::complex<double>(pars(180), pars(181));
   UV(4,2) = std::complex<double>(pars(182), pars(183));
   UV(4,3) = std::complex<double>(pars(184), pars(185));
   UV(4,4) = std::complex<double>(pars(186), pars(187));
   UV(4,5) = std::complex<double>(pars(188), pars(189));
   UV(5,0) = std::complex<double>(pars(190), pars(191));
   UV(5,1) = std::complex<double>(pars(192), pars(193));
   UV(5,2) = std::complex<double>(pars(194), pars(195));
   UV(5,3) = std::complex<double>(pars(196), pars(197));
   UV(5,4) = std::complex<double>(pars(198), pars(199));
   UV(5,5) = std::complex<double>(pars(200), pars(201));
   ZZ(0,0) = pars(202);
   ZZ(0,1) = pars(203);
   ZZ(1,0) = pars(204);
   ZZ(1,1) = pars(205);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(206);

   pars(22) = Re(Vd(0,0));
   pars(23) = Im(Vd(0,0));
   pars(24) = Re(Vd(0,1));
   pars(25) = Im(Vd(0,1));
   pars(26) = Re(Vd(0,2));
   pars(27) = Im(Vd(0,2));
   pars(28) = Re(Vd(1,0));
   pars(29) = Im(Vd(1,0));
   pars(30) = Re(Vd(1,1));
   pars(31) = Im(Vd(1,1));
   pars(32) = Re(Vd(1,2));
   pars(33) = Im(Vd(1,2));
   pars(34) = Re(Vd(2,0));
   pars(35) = Im(Vd(2,0));
   pars(36) = Re(Vd(2,1));
   pars(37) = Im(Vd(2,1));
   pars(38) = Re(Vd(2,2));
   pars(39) = Im(Vd(2,2));
   pars(40) = Re(Ud(0,0));
   pars(41) = Im(Ud(0,0));
   pars(42) = Re(Ud(0,1));
   pars(43) = Im(Ud(0,1));
   pars(44) = Re(Ud(0,2));
   pars(45) = Im(Ud(0,2));
   pars(46) = Re(Ud(1,0));
   pars(47) = Im(Ud(1,0));
   pars(48) = Re(Ud(1,1));
   pars(49) = Im(Ud(1,1));
   pars(50) = Re(Ud(1,2));
   pars(51) = Im(Ud(1,2));
   pars(52) = Re(Ud(2,0));
   pars(53) = Im(Ud(2,0));
   pars(54) = Re(Ud(2,1));
   pars(55) = Im(Ud(2,1));
   pars(56) = Re(Ud(2,2));
   pars(57) = Im(Ud(2,2));
   pars(58) = Re(Vu(0,0));
   pars(59) = Im(Vu(0,0));
   pars(60) = Re(Vu(0,1));
   pars(61) = Im(Vu(0,1));
   pars(62) = Re(Vu(0,2));
   pars(63) = Im(Vu(0,2));
   pars(64) = Re(Vu(1,0));
   pars(65) = Im(Vu(1,0));
   pars(66) = Re(Vu(1,1));
   pars(67) = Im(Vu(1,1));
   pars(68) = Re(Vu(1,2));
   pars(69) = Im(Vu(1,2));
   pars(70) = Re(Vu(2,0));
   pars(71) = Im(Vu(2,0));
   pars(72) = Re(Vu(2,1));
   pars(73) = Im(Vu(2,1));
   pars(74) = Re(Vu(2,2));
   pars(75) = Im(Vu(2,2));
   pars(76) = Re(Uu(0,0));
   pars(77) = Im(Uu(0,0));
   pars(78) = Re(Uu(0,1));
   pars(79) = Im(Uu(0,1));
   pars(80) = Re(Uu(0,2));
   pars(81) = Im(Uu(0,2));
   pars(82) = Re(Uu(1,0));
   pars(83) = Im(Uu(1,0));
   pars(84) = Re(Uu(1,1));
   pars(85) = Im(Uu(1,1));
   pars(86) = Re(Uu(1,2));
   pars(87) = Im(Uu(1,2));
   pars(88) = Re(Uu(2,0));
   pars(89) = Im(Uu(2,0));
   pars(90) = Re(Uu(2,1));
   pars(91) = Im(Uu(2,1));
   pars(92) = Re(Uu(2,2));
   pars(93) = Im(Uu(2,2));
   pars(94) = Re(Ve(0,0));
   pars(95) = Im(Ve(0,0));
   pars(96) = Re(Ve(0,1));
   pars(97) = Im(Ve(0,1));
   pars(98) = Re(Ve(0,2));
   pars(99) = Im(Ve(0,2));
   pars(100) = Re(Ve(1,0));
   pars(101) = Im(Ve(1,0));
   pars(102) = Re(Ve(1,1));
   pars(103) = Im(Ve(1,1));
   pars(104) = Re(Ve(1,2));
   pars(105) = Im(Ve(1,2));
   pars(106) = Re(Ve(2,0));
   pars(107) = Im(Ve(2,0));
   pars(108) = Re(Ve(2,1));
   pars(109) = Im(Ve(2,1));
   pars(110) = Re(Ve(2,2));
   pars(111) = Im(Ve(2,2));
   pars(112) = Re(Ue(0,0));
   pars(113) = Im(Ue(0,0));
   pars(114) = Re(Ue(0,1));
   pars(115) = Im(Ue(0,1));
   pars(116) = Re(Ue(0,2));
   pars(117) = Im(Ue(0,2));
   pars(118) = Re(Ue(1,0));
   pars(119) = Im(Ue(1,0));
   pars(120) = Re(Ue(1,1));
   pars(121) = Im(Ue(1,1));
   pars(122) = Re(Ue(1,2));
   pars(123) = Im(Ue(1,2));
   pars(124) = Re(Ue(2,0));
   pars(125) = Im(Ue(2,0));
   pars(126) = Re(Ue(2,1));
   pars(127) = Im(Ue(2,1));
   pars(128) = Re(Ue(2,2));
   pars(129) = Im(Ue(2,2));
   pars(130) = Re(UV(0,0));
   pars(131) = Im(UV(0,0));
   pars(132) = Re(UV(0,1));
   pars(133) = Im(UV(0,1));
   pars(134) = Re(UV(0,2));
   pars(135) = Im(UV(0,2));
   pars(136) = Re(UV(0,3));
   pars(137) = Im(UV(0,3));
   pars(138) = Re(UV(0,4));
   pars(139) = Im(UV(0,4));
   pars(140) = Re(UV(0,5));
   pars(141) = Im(UV(0,5));
   pars(142) = Re(UV(1,0));
   pars(143) = Im(UV(1,0));
   pars(144) = Re(UV(1,1));
   pars(145) = Im(UV(1,1));
   pars(146) = Re(UV(1,2));
   pars(147) = Im(UV(1,2));
   pars(148) = Re(UV(1,3));
   pars(149) = Im(UV(1,3));
   pars(150) = Re(UV(1,4));
   pars(151) = Im(UV(1,4));
   pars(152) = Re(UV(1,5));
   pars(153) = Im(UV(1,5));
   pars(154) = Re(UV(2,0));
   pars(155) = Im(UV(2,0));
   pars(156) = Re(UV(2,1));
   pars(157) = Im(UV(2,1));
   pars(158) = Re(UV(2,2));
   pars(159) = Im(UV(2,2));
   pars(160) = Re(UV(2,3));
   pars(161) = Im(UV(2,3));
   pars(162) = Re(UV(2,4));
   pars(163) = Im(UV(2,4));
   pars(164) = Re(UV(2,5));
   pars(165) = Im(UV(2,5));
   pars(166) = Re(UV(3,0));
   pars(167) = Im(UV(3,0));
   pars(168) = Re(UV(3,1));
   pars(169) = Im(UV(3,1));
   pars(170) = Re(UV(3,2));
   pars(171) = Im(UV(3,2));
   pars(172) = Re(UV(3,3));
   pars(173) = Im(UV(3,3));
   pars(174) = Re(UV(3,4));
   pars(175) = Im(UV(3,4));
   pars(176) = Re(UV(3,5));
   pars(177) = Im(UV(3,5));
   pars(178) = Re(UV(4,0));
   pars(179) = Im(UV(4,0));
   pars(180) = Re(UV(4,1));
   pars(181) = Im(UV(4,1));
   pars(182) = Re(UV(4,2));
   pars(183) = Im(UV(4,2));
   pars(184) = Re(UV(4,3));
   pars(185) = Im(UV(4,3));
   pars(186) = Re(UV(4,4));
   pars(187) = Im(UV(4,4));
   pars(188) = Re(UV(4,5));
   pars(189) = Im(UV(4,5));
   pars(190) = Re(UV(5,0));
   pars(191) = Im(UV(5,0));
   pars(192) = Re(UV(5,1));
   pars(193) = Im(UV(5,1));
   pars(194) = Re(UV(5,2));
   pars(195) = Im(UV(5,2));
   pars(196) = Re(UV(5,3));
   pars(197) = Im(UV(5,3));
   pars(198) = Re(UV(5,4));
   pars(199) = Im(UV(5,4));
   pars(200) = Re(UV(5,5));
   pars(201) = Im(UV(5,5));
   pars(202) = ZZ(0,0);
   pars(203) = ZZ(0,1);
   pars(204) = ZZ(1,0);
   pars(205) = ZZ(1,1);


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
   return "cSMHdCKMRHN";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   cSMHdCKMRHN_soft_parameters::run_to(scale, eps);
}







double CLASSNAME::get_mass_matrix_VG() const
{

   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{

   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

double CLASSNAME::get_mass_matrix_Hm() const
{

   const double mass_matrix_Hm = Re(0.25*(4*mu2 + 2*Lambdax*Sqr(v) + Sqr(g2)*
      Sqr(v)));

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{

   const auto mass_matrix_Hm = get_mass_matrix_Hm();
   MHm = mass_matrix_Hm;

   if (MHm < 0.) {
      problems.flag_running_tachyon(cSMHdCKMRHN_info::Hm);
   }

   MHm = AbsSqrt(MHm);
}

double CLASSNAME::get_mass_matrix_Ah() const
{

   const double mass_matrix_Ah = Re(0.25*(2*(2*mu2 + Lambdax*Sqr(v)) + Sqr(v)*
      Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{

   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_running_tachyon(cSMHdCKMRHN_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double CLASSNAME::get_mass_matrix_hh() const
{

   const double mass_matrix_hh = Re(mu2 + 1.5*Lambdax*Sqr(v));

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{

   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_running_tachyon(cSMHdCKMRHN_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(cSMHdCKMRHN_info::Fd, eigenvalue_error > precision *
      Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(cSMHdCKMRHN_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<std::complex<double>,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(cSMHdCKMRHN_info::Fe, eigenvalue_error > precision *
      Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<std::complex<double>,6,6> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(0,3) = 0.7071067811865475*v*Yv(0,0);
   mass_matrix_Fv(0,4) = 0.7071067811865475*v*Yv(1,0);
   mass_matrix_Fv(0,5) = 0.7071067811865475*v*Yv(2,0);
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(1,3) = 0.7071067811865475*v*Yv(0,1);
   mass_matrix_Fv(1,4) = 0.7071067811865475*v*Yv(1,1);
   mass_matrix_Fv(1,5) = 0.7071067811865475*v*Yv(2,1);
   mass_matrix_Fv(2,2) = 0;
   mass_matrix_Fv(2,3) = 0.7071067811865475*v*Yv(0,2);
   mass_matrix_Fv(2,4) = 0.7071067811865475*v*Yv(1,2);
   mass_matrix_Fv(2,5) = 0.7071067811865475*v*Yv(2,2);
   mass_matrix_Fv(3,3) = Mv(0,0);
   mass_matrix_Fv(3,4) = 0.5*Mv(0,1) + 0.5*Mv(1,0);
   mass_matrix_Fv(3,5) = 0.5*Mv(0,2) + 0.5*Mv(2,0);
   mass_matrix_Fv(4,4) = Mv(1,1);
   mass_matrix_Fv(4,5) = 0.5*Mv(1,2) + 0.5*Mv(2,1);
   mass_matrix_Fv(5,5) = Mv(2,2);

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   const auto mass_matrix_Fv(get_mass_matrix_Fv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV, eigenvalue_error);
   problems.flag_bad_mass(cSMHdCKMRHN_info::Fv, eigenvalue_error > precision *
      Abs(MFv(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Fv, MFv, UV);
#endif
   normalize_to_interval(UV);

}

double CLASSNAME::get_mass_matrix_VWp() const
{

   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void CLASSNAME::calculate_MVWp()
{

   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_running_tachyon(cSMHdCKMRHN_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void CLASSNAME::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> MVPVZ;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error);
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   
   double result = Re(mu2*v + 0.5*Cube(v)*Lambdax);

   return result;
}



double CLASSNAME::CphhHmconjHm() const
{
   
   const double result = -(v*Lambdax);

   return result;
}

double CLASSNAME::CpbargWpgZconjHm() const
{
   
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWpHm() const
{
   
   const double result = -0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpbargWpCgZHm() const
{
   
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWpCconjHm() const
{
   
   const double result = -0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpconjHmconjVWpVP() const
{
   
   const double result = -0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjHmconjVWpVZ() const
{
   
   const double result = 0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpAhAhHmconjHm() const
{
   
   const double result = -Lambdax;

   return result;
}

double CLASSNAME::CphhhhHmconjHm() const
{
   
   const double result = -Lambdax;

   return result;
}

double CLASSNAME::CpHmHmconjHmconjHm() const
{
   
   const double result = -2*Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhconjHmconjVWp() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double CLASSNAME::CphhconjHmconjVWp() const
{
   
   const double result = 0.5*g2;

   return result;
}

double CLASSNAME::CpHmconjHmVP() const
{
   
   const double result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpHmconjHmVZ() const
{
   
   const double result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpHmconjHmconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHmconjHmVZVZ() const
{
   
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())
      *Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(
      gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,
      Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvFeconjHmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(
      gI2,j1))*UV(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpFvFeconjHmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,
      Conj(UV(gI1,3 + j1))*Yv(j1,j2)));

   return result;
}

double CLASSNAME::CpAhAhhh() const
{
   
   const double result = -(v*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpbargWpgWpAh() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWpCgWpCAh() const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpAhAhAhAh() const
{
   
   const double result = -3*Lambdax;

   return result;
}

double CLASSNAME::CpAhAhhhhh() const
{
   
   const double result = -Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW(
      )) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhHmVWp() const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double CLASSNAME::CpAhAhconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ() const
{
   
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvAhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*(SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gI2,3 + j1))*UV(gI1,j2)) + SUM(
      j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(gI1,3 + j1))*UV(gI2,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvAhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *(SUM(j2,0,2,Conj(UV(gI2,j2))*SUM(j1,0,2,Conj(UV(gI1,3 + j1))*Yv(j1,j2))) +
      SUM(j2,0,2,Conj(UV(gI1,j2))*SUM(j1,0,2,Conj(UV(gI2,3 + j1))*Yv(j1,j2))));

   return result;
}

double CLASSNAME::Cphhhhhh() const
{
   
   const double result = -3*v*Lambdax;

   return result;
}

double CLASSNAME::CphhVZVZ() const
{
   
   const double result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CphhconjVWpVWp() const
{
   
   const double result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargWpgWphh() const
{
   
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargWpCgWpChh() const
{
   
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpbargZgZhh() const
{
   
   const double result = -0.25*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::Cphhhhhhhh() const
{
   
   const double result = -3*Lambdax;

   return result;
}

double CLASSNAME::CphhHmVWp() const
{
   
   const double result = -0.5*g2;

   return result;
}

double CLASSNAME::CphhhhconjVWpVWp() const
{
   
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ() const
{
   
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gI2,
      j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gI2,
      j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vu(gI2,
      j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvhhPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*(SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yv(j1,j2))*UV(gI2,3 + j1))*UV(gI1,j2)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(
      j1,j2))*UV(gI1,3 + j1))*UV(gI2,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvhhPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*(SUM(j2,0,2,Conj(UV(gI2
      ,j2))*SUM(j1,0,2,Conj(UV(gI1,3 + j1))*Yv(j1,j2))) + SUM(j2,0,2,Conj(UV(gI1,
      j2))*SUM(j1,0,2,Conj(UV(gI2,3 + j1))*Yv(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbargGgGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpbarFdFdVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFdFdVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   
   const double result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   
   const double result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpHmVPVWp() const
{
   
   const double result = -0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpgWpVP() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgWpCVP() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpHmconjHmVPVP() const
{
   
   const std::complex<double> result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834*g1
      *Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpconjVWpVPVWp() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPR(int gI1, int gI2) const
{
   
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFeFeVPPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVPPR(int gI1, int gI2) const
{
   
   const double result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFuFuVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVPPR(int gI1, int gI2) const
{
   
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpconjVWpVPVPVWp3() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVPVPVWp1() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVPVPVWp2() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHmVWpVZ() const
{
   
   const double result = 0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpgWpVZ() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgWpCVZ() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZ() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbarFdFdVZPL(int gI1, int gI2) const
{
   
   const double result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPR(int gI1, int gI2) const
{
   
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFeFeVZPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) -
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPR(int gI1, int gI2) const
{
   
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW(
      ));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(UV(gI2,j1))*UV(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpFvFvVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*
      g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(UV(gI1,j1))*UV(gI2,j1));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZVZ1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZVZ2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpVWpVZVZ3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargPgWpconjVWp() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWpCgPconjVWp() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWpCgZconjVWp() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpbargZgWpconjVWp() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuconjVWpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vu(
      gI2,j1))*Vd(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFdFuconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvconjVWpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(UV(
      gI2,j1))*Ve(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFeFvconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp2() const
{
   
   const double result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp1() const
{
   
   const double result = 2*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpconjVWpVWpVWp3() const
{
   
   const double result = -Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(gI1,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(Vd(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yd(j1,gO1))*Ud(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Ud(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vd(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(ThetaW
      ())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(Vd(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Sin(ThetaW(
      )),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Sin(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(Vu(gI2,j2))*Yd(
      gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(
      gI2,j1)),0);

   return result;
}

double CLASSNAME::CpbarUFdFuconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuconjVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vu(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(gI1,j1)),0);

   return result;
}

double CLASSNAME::CpbarUFuFdVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Vd(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vd(gI2,j2))*Yu(
      gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(
      gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(Vu(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yu(j1,gO1))*Uu(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Uu(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vu(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5163977794943222*g1*Cos(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(Vu(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(ThetaW
      ())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI1,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(Ve(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Ye(j1,gO1))*Ue(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(ThetaW
      ())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(Ve(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW(
      )),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Sin(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(UV(gI2,j2))*Ye(
      gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yv(j1,gO1))*UV(
      gI2,3 + j1)),0);

   return result;
}

double CLASSNAME::CpbarUFeFvconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvconjVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(UV(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUFvHmPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(
      j1,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUFvHmPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*
      KroneckerDelta(gO1,3 + j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvAhPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(UV(gI1,3 + j1))*Yv(j1,gO2)),0) + std::
      complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(UV(gI1,j2))*SUM(j1,0,
      2,KroneckerDelta(gO2,3 + j1)*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvAhPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j1,0,2,Conj(Yv(j1,gO1))*UV(gI1,3 + j1)),0) - std::
      complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))
      *KroneckerDelta(gO1,3 + j1))*UV(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUFvconjVWpPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*g2*Ve(gI2,
      gO2),0);

   return result;
}

double CLASSNAME::CpbarFeUFvconjVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpUFvFeVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUFvFeVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(Ve(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpUFvFeconjHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUFvFeconjHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(
      gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvhhPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(UV(gI2,3 + j1))*Yv(j1,gO2)),0) - 0.7071067811865475*SUM(j2,0,2,Conj(UV(
      gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvhhPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yv(j1,gO1))*UV(gI2,3 + j1)),0) - 0.7071067811865475*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yv(j1,j2))*KroneckerDelta(gO1,3 + j1))*UV(gI2,j2));

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Cos(ThetaW())*UV(gI2,gO2
      ),0) + IF(gO2 < 3,-0.3872983346207417*g1*Sin(ThetaW())*UV(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpFvUFvVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5*g2*Conj(UV(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gO1 < 3,0.3872983346207417*g1*Conj(UV(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj
      (Ud(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(
      gI2,j1))*Vd(gO1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHmPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = SUM(j2,0,2,Conj(UV(gI2,j2))*SUM(j1,0,2,Conj
      (Ue(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHmPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*UV(
      gI2,3 + j1))*Ve(gO1,j2));

   return result;
}

double CLASSNAME::CpbarFuFdVWpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdVWpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vd(
      gI2,j1))*Vu(gO1,j1));

   return result;
}


std::complex<double> CLASSNAME::self_energy_Hm_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhHmconjHm();
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVZ))*CpbargWpCgZHm()*CpbargZgWpCconjHm());
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVWp))*CpbargWpgZconjHm()*CpbargZgWpHm());
   result += 2*AbsSqr(CpconjHmconjVWpVP())*(-1 + 2*B0(Sqr(p),0,Sqr(MVWp)));
   result += 2*AbsSqr(CpconjHmconjVWpVZ())*(-1 + 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVWp)));
   result += -0.5*A0(Sqr(Mhh))*CphhhhHmconjHm();
   result += AbsSqr(CphhHmconjHm())*B0(Sqr(p),Sqr(MHm),Sqr(Mhh));
   result += -(A0(Sqr(MHm))*CpHmHmconjHmconjHm());
   result += AbsSqr(CpAhconjHmconjVWp())*F0(Sqr(p),Sqr(MAh),Sqr(MVWp));
   result += AbsSqr(CphhconjHmconjVWp())*F0(Sqr(p),Sqr(Mhh),Sqr(MVWp));
   result += AbsSqr(CpHmconjHmVP())*F0(Sqr(p),Sqr(MHm),0);
   result += AbsSqr(CpHmconjHmVZ())*F0(Sqr(p),Sqr(MHm),Sqr(MVZ));
   result += 4*A0(Sqr(MVWp))*CpHmconjHmconjVWpVWp() - 2*CpHmconjHmconjVWpVWp()*Sqr
      (MVWp);
   result += CpHmconjHmVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFdconjHmPL(gI1,gI2)) +
      AbsSqr(CpbarFuFdconjHmPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))))
      ;
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFuFdconjHmPR(gI1,gI2))*CpbarFuFdconjHmPL(gI1,gI2) + Conj(
      CpbarFuFdconjHmPL(gI1,gI2))*CpbarFuFdconjHmPR(gI1,gI2))*MFd(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,(AbsSqr(CpFvFeconjHmPL(gI1,gI2)) + AbsSqr(
      CpFvFeconjHmPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))));
   result += -2*SUM(gI1,0,5,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpFvFeconjHmPR(gI1,gI2))*CpFvFeconjHmPL(gI1,gI2) + Conj(
      CpFvFeconjHmPL(gI1,gI2))*CpFvFeconjHmPR(gI1,gI2))*MFe(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Ah_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhAhAh();
   result += AbsSqr(CpAhAhhh())*B0(Sqr(p),Sqr(Mhh),Sqr(MAh));
   result += -0.5*A0(Sqr(Mhh))*CpAhAhhhhh();
   result += -(A0(Sqr(MHm))*CpAhAhHmconjHm());
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpCgWpCAh()));
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpgWpAh()));
   result += AbsSqr(CpAhhhVZ())*F0(Sqr(p),Sqr(Mhh),Sqr(MVZ));
   result += 2*AbsSqr(CpAhHmVWp())*F0(Sqr(p),Sqr(MHm),Sqr(MVWp));
   result += 4*A0(Sqr(MVWp))*CpAhAhconjVWpVWp() - 2*CpAhAhconjVWpVWp()*Sqr(MVWp);
   result += CpAhAhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdAhPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeAhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuAhPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdAhPR(gI1,gI2))*CpbarFdFdAhPL(gI1,gI2) + Conj(
      CpbarFdFdAhPL(gI1,gI2))*CpbarFdFdAhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeAhPR(gI1,gI2))*CpbarFeFeAhPL(gI1,gI2) + Conj(
      CpbarFeFeAhPL(gI1,gI2))*CpbarFeFeAhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuAhPR(gI1,gI2))*CpbarFuFuAhPL(gI1,gI2) + Conj(
      CpbarFuFuAhPL(gI1,gI2))*CpbarFuFuAhPR(gI1,gI2))*MFu(gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(AbsSqr(CpFvFvAhPL(gI1,gI2)) + AbsSqr(
      CpFvFvAhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2)))));
   result += -SUM(gI1,0,5,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2
      )))*(Conj(CpFvFvAhPR(gI1,gI2))*CpFvFvAhPL(gI1,gI2) + Conj(CpFvFvAhPL(gI1,gI2
      ))*CpFvFvAhPR(gI1,gI2))*MFv(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_hh_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CpAhAhhh())*B0(Sqr(p),Sqr(MAh),Sqr(MAh));
   result += -0.5*A0(Sqr(MAh))*CpAhAhhhhh();
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpCgWpChh()));
   result += -(B0(Sqr(p),Sqr(MVWp),Sqr(MVWp))*Sqr(CpbargWpgWphh()));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*Sqr(CpbargZgZhh()));
   result += 2*AbsSqr(CphhconjVWpVWp())*(-1 + 2*B0(Sqr(p),Sqr(MVWp),Sqr(MVWp)));
   result += 0.5*AbsSqr(Cphhhhhh())*B0(Sqr(p),Sqr(Mhh),Sqr(Mhh));
   result += -0.5*A0(Sqr(Mhh))*Cphhhhhhhh();
   result += -(A0(Sqr(MHm))*CphhhhHmconjHm());
   result += AbsSqr(CphhHmconjHm())*B0(Sqr(p),Sqr(MHm),Sqr(MHm));
   result += AbsSqr(CphhVZVZ())*(-1 + 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ)));
   result += AbsSqr(CpAhhhVZ())*F0(Sqr(p),Sqr(MAh),Sqr(MVZ));
   result += 2*AbsSqr(CphhHmVWp())*F0(Sqr(p),Sqr(MHm),Sqr(MVWp));
   result += 4*A0(Sqr(MVWp))*CphhhhconjVWpVWp() - 2*CphhhhconjVWpVWp()*Sqr(MVWp);
   result += CphhhhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdhhPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdhhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFehhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFehhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuhhPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuhhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdhhPR(gI1,gI2))*CpbarFdFdhhPL(gI1,gI2) + Conj(
      CpbarFdFdhhPL(gI1,gI2))*CpbarFdFdhhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFehhPR(gI1,gI2))*CpbarFeFehhPL(gI1,gI2) + Conj(
      CpbarFeFehhPL(gI1,gI2))*CpbarFeFehhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuhhPR(gI1,gI2))*CpbarFuFuhhPL(gI1,gI2) + Conj(
      CpbarFuFuhhPL(gI1,gI2))*CpbarFuFuhhPR(gI1,gI2))*MFu(gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(AbsSqr(CpFvFvhhPL(gI1,gI2)) + AbsSqr(
      CpFvFvhhPR(gI1,gI2)))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2)))));
   result += -SUM(gI1,0,5,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2
      )))*(Conj(CpFvFvhhPR(gI1,gI2))*CpFvFvhhPL(gI1,gI2) + Conj(CpFvFvhhPL(gI1,gI2
      ))*CpFvFvhhPR(gI1,gI2))*MFv(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VG_1loop(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpbargGgGVG())*B00(Sqr(p),Sqr(MVG),Sqr(MVG));
   result += -(AbsSqr(CpVGVGVG())*(15*B00(Sqr(p),0,0) + Sqr(p) + 6*B0(Sqr(p),0,0)*
      Sqr(p)));
   result += 0;
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVGPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVGPL(gI1,
      gI2))*CpbarFdFdVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVGPL(gI1,
      gI2))*CpbarFuFuVGPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VP_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWpCgWpCVP())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += AbsSqr(CpbargWpgWpVP())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += -(A0(Sqr(MVWp))*(CpconjVWpVPVPVWp1() + CpconjVWpVPVPVWp2() + 4*
      CpconjVWpVPVPVWp3()));
   result += -4*AbsSqr(CpHmconjHmVP())*B00(Sqr(p),Sqr(MHm),Sqr(MHm));
   result += A0(Sqr(MHm))*CpHmconjHmVPVP();
   result += 2*AbsSqr(CpHmVPVWp())*B0(Sqr(p),Sqr(MVWp),Sqr(MHm));
   result += 2*CpconjVWpVPVPVWp3()*Sqr(MVWp);
   result += -0.6666666666666666*AbsSqr(CpconjVWpVPVWp())*(3*A0(Sqr(MVWp)) + 15*
      B00(Sqr(p),Sqr(MVWp),Sqr(MVWp)) - 6*Sqr(MVWp) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWp),Sqr(MVWp))*(Sqr(MVWp) + 2*Sqr(p)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,
      gI2))*CpbarFdFdVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,
      gI2))*CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,
      gI2))*CpbarFuFuVPPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VZ_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*A0(Sqr(MAh))*CpAhAhVZVZ();
   result += -4*AbsSqr(CpAhhhVZ())*B00(Sqr(p),Sqr(MAh),Sqr(Mhh));
   result += AbsSqr(CpbargWpCgWpCVZ())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += AbsSqr(CpbargWpgWpVZ())*B00(Sqr(p),Sqr(MVWp),Sqr(MVWp));
   result += -(A0(Sqr(MVWp))*(4*CpconjVWpVWpVZVZ1() + CpconjVWpVWpVZVZ2() +
      CpconjVWpVWpVZVZ3()));
   result += 0.5*A0(Sqr(Mhh))*CphhhhVZVZ();
   result += AbsSqr(CphhVZVZ())*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh));
   result += -4*AbsSqr(CpHmconjHmVZ())*B00(Sqr(p),Sqr(MHm),Sqr(MHm));
   result += A0(Sqr(MHm))*CpHmconjHmVZVZ();
   result += 2*AbsSqr(CpHmVWpVZ())*B0(Sqr(p),Sqr(MVWp),Sqr(MHm));
   result += 2*CpconjVWpVWpVZVZ1()*Sqr(MVWp);
   result += -0.6666666666666666*AbsSqr(CpconjVWpVWpVZ())*(3*A0(Sqr(MVWp)) + 15*
      B00(Sqr(p),Sqr(MVWp),Sqr(MVWp)) - 6*Sqr(MVWp) + Sqr(p) + 3*B0(Sqr(p),Sqr(
      MVWp),Sqr(MVWp))*(Sqr(MVWp) + 2*Sqr(p)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,
      gI2))*CpbarFdFdVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZPL(gI1,
      gI2))*CpbarFeFeVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZPL(gI1,
      gI2))*CpbarFuFuVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(AbsSqr(CpFvFvVZPL(gI1,gI2)) + AbsSqr(
      CpFvFvVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(Sqr(p),
      Sqr(MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpFvFvVZPL(gI1,gI2))*
      CpFvFvVZPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VWp_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*A0(Sqr(MAh))*CpAhAhconjVWpVWp();
   result += -4*AbsSqr(CpAhconjHmconjVWp())*B00(Sqr(p),Sqr(MAh),Sqr(MHm));
   result += AbsSqr(CpbargPgWpconjVWp())*B00(Sqr(p),Sqr(MVWp),Sqr(MVP));
   result += AbsSqr(CpbargWpCgPconjVWp())*B00(Sqr(p),Sqr(MVP),Sqr(MVWp));
   result += AbsSqr(CpbargWpCgZconjVWp())*B00(Sqr(p),Sqr(MVZ),Sqr(MVWp));
   result += AbsSqr(CpbargZgWpconjVWp())*B00(Sqr(p),Sqr(MVWp),Sqr(MVZ));
   result += AbsSqr(CpconjHmconjVWpVP())*B0(Sqr(p),0,Sqr(MHm));
   result += AbsSqr(CpconjHmconjVWpVZ())*B0(Sqr(p),Sqr(MVZ),Sqr(MHm));
   result += -(A0(Sqr(MVWp))*(CpconjVWpconjVWpVWpVWp1() + 4*
      CpconjVWpconjVWpVWpVWp2() + CpconjVWpconjVWpVWpVWp3()));
   result += 0;
   result += -4*AbsSqr(CphhconjHmconjVWp())*B00(Sqr(p),Sqr(Mhh),Sqr(MHm));
   result += AbsSqr(CphhconjVWpVWp())*B0(Sqr(p),Sqr(MVWp),Sqr(Mhh));
   result += 0.5*A0(Sqr(Mhh))*CphhhhconjVWpVWp();
   result += A0(Sqr(MHm))*CpHmconjHmconjVWpVWp();
   result += 2*CpconjVWpconjVWpVWpVWp2()*Sqr(MVWp);
   result += -(AbsSqr(CpconjVWpVPVWp())*(A0(Sqr(MVWp)) + 10*B00(Sqr(p),Sqr(MVWp),0
      ) - 2*Sqr(MVWp) + 0.6666666666666666*Sqr(p) + B0(Sqr(p),Sqr(MVWp),0)*(Sqr(
      MVWp) + 4*Sqr(p))));
   result += -0.5*A0(Sqr(MVZ))*(4*CpconjVWpVWpVZVZ1() + CpconjVWpVWpVZVZ2() +
      CpconjVWpVWpVZVZ3()) + CpconjVWpVWpVZVZ1()*Sqr(MVZ);
   result += -(AbsSqr(CpconjVWpVWpVZ())*(A0(Sqr(MVWp)) + A0(Sqr(MVZ)) + 10*B00(Sqr
      (p),Sqr(MVZ),Sqr(MVWp)) - 2*(Sqr(MVWp) + Sqr(MVZ) - 0.3333333333333333*Sqr(p
      )) + B0(Sqr(p),Sqr(MVZ),Sqr(MVWp))*(Sqr(MVWp) + Sqr(MVZ) + 4*Sqr(p))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFuconjVWpPL(gI1,gI2)) +
      AbsSqr(CpbarFdFuconjVWpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(gI2)))
      + 4*B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFu(gI2)))*MFd(gI1)*MFu(gI2)*Re(Conj(
      CpbarFdFuconjVWpPL(gI1,gI2))*CpbarFdFuconjVWpPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,5,(AbsSqr(CpbarFeFvconjVWpPL(gI1,gI2)) + AbsSqr
      (CpbarFeFvconjVWpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(gI2))) + 4*B0
      (Sqr(p),Sqr(MFe(gI1)),Sqr(MFv(gI2)))*MFe(gI1)*MFv(gI2)*Re(Conj(
      CpbarFeFvconjVWpPL(gI1,gI2))*CpbarFeFvconjVWpPR(gI1,gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarUFdFdAhPL(gO2
      ,gI1))*CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarUFdFdhhPL(gO2
      ,gI2))*CpbarUFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),0))*Conj(
      CpbarUFdFdVPPR(gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2)*MFu(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(CpbarUFdFuHmPL(gO2
      ,gI2))*CpbarUFdFuHmPR(gO1,gI2)*MFu(gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(
      CpbarUFdFdAhPR(gO2,gI1))*CpbarUFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(
      CpbarUFdFdhhPR(gO2,gI2))*CpbarUFdFdhhPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPL(
      gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPL(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPL(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(
      CpbarUFdFuHmPR(gO2,gI2))*CpbarUFdFuHmPR(gO1,gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(
      CpbarUFdFdAhPL(gO2,gI1))*CpbarUFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(
      CpbarUFdFdhhPL(gO2,gI2))*CpbarUFdFdhhPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPR(
      gO2,gI2))*CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(
      CpbarUFdFuHmPL(gO2,gI2))*CpbarUFdFuHmPL(gO1,gI2));

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
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(CpbarUFuFdconjHmPL
      (gO2,gI2))*CpbarUFuFdconjHmPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(
      CpbarUFuFuVPPR(gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPR(gO2,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarUFuFdconjHmPR(gO2,gI2))*CpbarUFuFdconjHmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPR(gO2,gI2))*CpbarUFuFuhhPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPL(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPL(gO2,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarUFuFdconjHmPL(gO2,gI2))*CpbarUFuFdconjHmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPL(gO2,gI2))*CpbarUFuFuhhPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

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
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarUFeFeAhPL(gO2
      ,gI1))*CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarUFeFehhPL(gO2
      ,gI2))*CpbarUFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),0))*Conj(
      CpbarUFeFeVPPR(gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,5,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2)*MFv(gI2));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(CpbarUFeFvHmPL(gO2
      ,gI2))*CpbarUFeFvHmPR(gO1,gI2)*MFv(gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(
      CpbarUFeFeAhPR(gO2,gI1))*CpbarUFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(
      CpbarUFeFehhPR(gO2,gI2))*CpbarUFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPL(
      gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPL(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPL(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(
      CpbarUFeFvHmPR(gO2,gI2))*CpbarUFeFvHmPR(gO1,gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(
      CpbarUFeFeAhPL(gO2,gI1))*CpbarUFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(
      CpbarUFeFehhPL(gO2,gI2))*CpbarUFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPR(
      gO2,gI2))*CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(
      CpbarUFeFvHmPL(gO2,gI2))*CpbarUFeFvHmPL(gO1,gI2));

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
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MHm))*Conj(CpbarFeUFvHmPL(gI1
      ,gO2))*CpbarFeUFvHmPR(gI1,gO1)*MFe(gI1));
   result += SUM(gI1,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh))*Conj(CpFvUFvAhPL(gI1,
      gO2))*CpFvUFvAhPR(gI1,gO1)*MFv(gI1));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeUFvconjVWpPL(gI2,gO2))*CpbarFeUFvconjVWpPR(gI2,gO1)*MFe(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm))*Conj(CpUFvFeconjHmPL(
      gO2,gI2))*CpUFvFeconjHmPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpUFvFeVWpPR(gO2,gI2))*CpUFvFeVWpPL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh))*Conj(CpFvUFvhhPL(gI2,
      gO2))*CpFvUFvhhPR(gI2,gO1)*MFv(gI2));
   result += -4*SUM(gI2,0,5,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpFvUFvVZPL(gI2,gO2))*CpFvUFvVZPR(gI2,gO1)*MFv(gI2));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MHm))*Conj(
      CpbarFeUFvHmPR(gI1,gO2))*CpbarFeUFvHmPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh))*Conj(CpFvUFvAhPR(
      gI1,gO2))*CpFvUFvAhPR(gI1,gO1));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeUFvconjVWpPR(gI2,gO2))*CpbarFeUFvconjVWpPR(gI2,gO1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm))*Conj(
      CpUFvFeconjHmPR(gO2,gI2))*CpUFvFeconjHmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpUFvFeVWpPL(gO2,gI2))*CpUFvFeVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh))*Conj(CpFvUFvhhPR(
      gI2,gO2))*CpFvUFvhhPR(gI2,gO1));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpFvUFvVZPR(gI2,gO2))*CpFvUFvVZPR(gI2,gO1));

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
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MHm))*Conj(
      CpbarFeUFvHmPL(gI1,gO2))*CpbarFeUFvHmPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh))*Conj(CpFvUFvAhPL(
      gI1,gO2))*CpFvUFvAhPL(gI1,gO1));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeUFvconjVWpPL(gI2,gO2))*CpbarFeUFvconjVWpPL(gI2,gO1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHm))*Conj(
      CpUFvFeconjHmPL(gO2,gI2))*CpUFvFeconjHmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWp)))*Conj(
      CpUFvFeVWpPR(gO2,gI2))*CpUFvFeVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh))*Conj(CpFvUFvhhPL(
      gI2,gO2))*CpFvUFvhhPL(gI2,gO1));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ)))*Conj(
      CpFvUFvVZPL(gI2,gO2))*CpFvUFvVZPL(gI2,gO1));

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

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPL(gO2,
      gI1))*CpbarFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPL(gO2,
      gI2))*CpbarFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2)*MFu(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(CpbarFdFuHmPL(gO2,
      gI2))*CpbarFdFuHmPR(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPR
      (gO2,gI1))*CpbarFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPR
      (gO2,gI2))*CpbarFdFdhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPL(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPL(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(CpbarFdFuHmPR
      (gO2,gI2))*CpbarFdFuHmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh))*Conj(CpbarFdFdAhPL
      (gO2,gI1))*CpbarFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh))*Conj(CpbarFdFdhhPL
      (gO2,gI2))*CpbarFdFdhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ)))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWp)))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHm))*Conj(CpbarFdFuHmPL
      (gO2,gI2))*CpbarFdFuHmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPL(gO2,
      gI1))*CpbarFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPL(gO2,
      gI2))*CpbarFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,5,(-0.5 + B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2)*MFv(gI2));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(CpbarFeFvHmPL(gO2,
      gI2))*CpbarFeFvHmPR(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPR
      (gO2,gI1))*CpbarFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPR
      (gO2,gI2))*CpbarFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPL(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPL(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(CpbarFeFvHmPR
      (gO2,gI2))*CpbarFeFvHmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh))*Conj(CpbarFeFeAhPL
      (gO2,gI1))*CpbarFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh))*Conj(CpbarFeFehhPL
      (gO2,gI2))*CpbarFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ)))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,5,(0.5 + B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWp)))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHm))*Conj(CpbarFeFvHmPL
      (gO2,gI2))*CpbarFeFvHmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPL(gO2,
      gI1))*CpbarFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(CpbarFuFdconjHmPL(
      gO2,gI2))*CpbarFuFdconjHmPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPL(gO2,
      gI2))*CpbarFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR
      (gO2,gI2))*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPR
      (gO2,gI1))*CpbarFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarFuFdconjHmPR(gO2,gI2))*CpbarFuFdconjHmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPL(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPR
      (gO2,gI2))*CpbarFuFuhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPL(
      gO2,gI2))*CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPL(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarFuFuAhPL
      (gO2,gI1))*CpbarFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarFuFdconjHmPL(gO2,gI2))*CpbarFuFdconjHmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarFuFuhhPL
      (gO2,gI2))*CpbarFuFuhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR(
      gO2,gI2))*CpbarFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(CpbarUFuFdconjHmPL
      (gO2,gI2))*CpbarUFuFdconjHmPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),0))*Conj(
      CpbarUFuFuVPPR(gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPR(gO2,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarUFuFdconjHmPR(gO2,gI2))*CpbarUFuFdconjHmPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPR(gO2,gI2))*CpbarUFuFuhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPL(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh))*Conj(
      CpbarUFuFuAhPL(gO2,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHm))*Conj(
      CpbarUFuFdconjHmPL(gO2,gI2))*CpbarUFuFdconjHmPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWp)))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh))*Conj(
      CpbarUFuFuhhPL(gO2,gI2))*CpbarUFuFuhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ)))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::tadpole_hh_1loop() const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MAh))*CpAhAhhh();
   result += A0(Sqr(MVWp))*CpbargWpCgWpChh();
   result += A0(Sqr(MVWp))*CpbargWpgWphh();
   result += A0(Sqr(MVZ))*CpbargZgZhh();
   result += -0.5*A0(Sqr(Mhh))*Cphhhhhh();
   result += -(A0(Sqr(MHm))*CphhHmconjHm());
   result += 4*A0(Sqr(MVWp))*CphhconjVWpVWp() - 2*CphhconjVWpVWp()*Sqr(MVWp);
   result += CphhVZVZ()*(2*A0(Sqr(MVZ)) - Sqr(MVZ));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdhhPL(gI1,gI1) +
      CpbarFdFdhhPR(gI1,gI1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFehhPL(gI1,gI1) +
      CpbarFeFehhPR(gI1,gI1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuhhPL(gI1,gI1) +
      CpbarFuFuhhPR(gI1,gI1))*MFu(gI1));
   result += SUM(gI1,0,5,A0(Sqr(MFv(gI1)))*(CpFvFvhhPL(gI1,gI1) + CpFvFvhhPR(gI1,
      gI1))*MFv(gI1));

   return result * oneOver16PiSqr;
}













void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(cSMHdCKMRHN_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const double M_tree(get_mass_matrix_hh());
      const double p = old_Mhh;
      double self_energy = Re(self_energy_hh_1loop(p));
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(Mhh) = SignedAbsSqrt(mass_sqr);

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(cSMHdCKMRHN_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(cSMHdCKMRHN_info::hh);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(cSMHdCKMRHN_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
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
      problems.flag_bad_mass(cSMHdCKMRHN_info::Fd, eigenvalue_error > precision
          * Abs(eigen_values(0)));
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
   double qcd_1l = 0.;

   {
      const double currentScale = get_scale();
      qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(2))/Sqr(currentScale)
         ))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = -0.005284774766427138*Quad(g3) - 0.0032348537833770956*Log(Sqr(
         currentScale)/Sqr(MFu(2)))*Quad(g3) - 0.0008822328500119351*Quad(g3)*
         Sqr(Log(Sqr(currentScale)/Sqr(MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = -0.00003352082872926087*Power6(g3)*(35.70257721711602 + 1.*Cube(
         Log(Sqr(currentScale)/Sqr(MFu(2)))) + 15.3874108148848*Log(Sqr(
         currentScale)/Sqr(MFu(2))) + 5.378787878787879*Sqr(Log(Sqr(
         currentScale)/Sqr(MFu(2)))));
   }

   const Eigen::Matrix<std::complex<double>,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      Eigen::Matrix<std::complex<double>,3,3> self_energy_1;
      Eigen::Matrix<std::complex<double>,3,3> self_energy_PL;
      Eigen::Matrix<std::complex<double>,3,3> self_energy_PR;
      for (int i1 = 0; i1 < 3; ++i1) {
         for (int i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = self_energy_Fu_1loop_1_heavy(p,i1,i2);
               self_energy_PL(i1,i2) = self_energy_Fu_1loop_PL_heavy(p,i1,i2);
               self_energy_PR(i1,i2) = self_energy_Fu_1loop_PR_heavy(p,i1,i2);
            } else {
               self_energy_1(i1,i2)  = self_energy_Fu_1loop_1(p,i1,i2);
               self_energy_PL(i1,i2) = self_energy_Fu_1loop_PL(p,i1,i2);
               self_energy_PR(i1,i2) = self_energy_Fu_1loop_PR(p,i1,i2);
            }
         }
      }
      Eigen::Matrix<std::complex<double>,3,3> delta_M(- self_energy_PR * M_tree
          - M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
      const Eigen::Matrix<std::complex<double>,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(cSMHdCKMRHN_info::Fu, eigenvalue_error > precision
          * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
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
      problems.flag_bad_mass(cSMHdCKMRHN_info::Fe, eigenvalue_error > precision
          * Abs(eigen_values(0)));
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

void CLASSNAME::calculate_MVWp_pole()
{
   if (!force_output && problems.is_running_tachyon(cSMHdCKMRHN_info::VWp))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWp));
   const double p = MVWp;
   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VWp);

   PHYSICAL(MVWp) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(cSMHdCKMRHN_info::VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = Sqr(MVWp) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VWp);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(cSMHdCKMRHN_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VZ);

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::calculate_MFv_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fv_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fv_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fv_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL +
      self_energy_PR);

   return m_drbar;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mf_1loop = - self_energy_1/m_sm_drbar - self_energy_PL -
      self_energy_PR;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0.;

   qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(idx))/Sqr(currentScale))
      )*Sqr(g3);

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 0.005284774766427138*Quad(g3) + 0.0032348537833770956
         *Log(Sqr(currentScale)/Sqr(MFu(idx)))*Quad(g3) + 0.0008822328500119351
         *Quad(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(idx))));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   if (get_thresholds() > 2 && threshold_corrections.mt > 2) {
      qcd_3l = -0.0008783313853540776*Power6(g3) - 5.078913443827405e-6*Cube(
         Log(Sqr(currentScale)/Sqr(MFu(idx))))*Power6(g3) -
         0.00011624292515063117*Log(Sqr(currentScale)/Sqr(MFu(idx)))*Power6(g3)
         - 0.00002183932780845784*Power6(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(
         idx))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL -
      self_energy_PR;
   double qcd_2l = 0.;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop + qcd_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(cSMHdCKMRHN_info::VWp);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}



std::ostream& operator<<(std::ostream& ostr, const cSMHdCKMRHN_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

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

// File generated at Mon 5 Nov 2018 12:47:51

/**
 * @file cSMHdCKM_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 5 Nov 2018 12:47:51 with FlexibleSUSY
 * 2.2.0 (git commit: fff0745460c710cd894e99b5b50967ac42dc9aba) and SARAH 4.14.0 .
 */

#ifndef cSMHdCKM_MASS_EIGENSTATES_H
#define cSMHdCKM_MASS_EIGENSTATES_H

#include "cSMHdCKM_info.hpp"
#include "cSMHdCKM_physical.hpp"
#include "cSMHdCKM_soft_parameters.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>
#include <string>

#include <Eigen/Core>

namespace flexiblesusy {

class cSMHdCKM_ewsb_solver_interface;
/**
 * @class cSMHdCKM_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class cSMHdCKM_mass_eigenstates : public cSMHdCKM_soft_parameters {
public:
   explicit cSMHdCKM_mass_eigenstates(const cSMHdCKM_input_parameters& input_ = cSMHdCKM_input_parameters());
   cSMHdCKM_mass_eigenstates(const cSMHdCKM_mass_eigenstates&) = default;
   cSMHdCKM_mass_eigenstates(cSMHdCKM_mass_eigenstates&&) = default;
   virtual ~cSMHdCKM_mass_eigenstates() = default;
   cSMHdCKM_mass_eigenstates& operator=(const cSMHdCKM_mass_eigenstates&) = default;
   cSMHdCKM_mass_eigenstates& operator=(cSMHdCKM_mass_eigenstates&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear() override;
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   Eigen::ArrayXd get_DRbar_masses_and_mixings() const;
   Eigen::ArrayXd get_extra_parameters() const;
   void do_calculate_sm_pole_masses(bool);
   bool do_calculate_sm_pole_masses() const;
   void do_calculate_bsm_pole_masses(bool);
   bool do_calculate_bsm_pole_masses() const;
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_DRbar_masses(const Eigen::ArrayXd&);
   void set_DRbar_masses_and_mixings(const Eigen::ArrayXd&);
   void set_extra_parameters(const Eigen::ArrayXd&);
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const cSMHdCKM_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const cSMHdCKM_physical& get_physical() const;
   cSMHdCKM_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   void set_ewsb_solver(const std::shared_ptr<cSMHdCKM_ewsb_solver_interface>&);
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   void calculate_spectrum();
   void clear_problems();
   std::string name() const;
   void run_to(double scale, double eps = -1.0) override;
   void print(std::ostream& out = std::cerr) const override;
   void set_precision(double);
   double get_precision() const;


   double get_MVG() const { return MVG; }
   double get_MHm() const { return MHm; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_MAh() const { return MAh; }
   double get_Mhh() const { return Mhh; }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_MVWp() const { return MVWp; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }

   

   
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   std::complex<double> get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   std::complex<double> get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   std::complex<double> get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   std::complex<double> get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   std::complex<double> get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   std::complex<double> get_Ue(int i, int k) const { return Ue(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const { return ZZ; }
   double get_ZZ(int i, int k) const { return ZZ(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_UV() const { return UV; }
   std::complex<double> get_UV(int i, int k) const { return UV(i,k); }



   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_Hm() const;
   void calculate_MHm();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   double get_mass_matrix_Ah() const;
   void calculate_MAh();
   double get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<std::complex<double>,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   double get_mass_matrix_VWp() const;
   void calculate_MVWp();
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const;
   void calculate_MVPVZ();

   double get_ewsb_eq_hh_1() const;

   double CphhHmconjHm() const;
   double CpbargWpgZconjHm() const;
   double CpbargZgWpHm() const;
   double CpbargWpCgZHm() const;
   double CpbargZgWpCconjHm() const;
   double CpconjHmconjVWpVP() const;
   double CpconjHmconjVWpVZ() const;
   double CpAhAhHmconjHm() const;
   double CphhhhHmconjHm() const;
   double CpHmHmconjHmconjHm() const;
   std::complex<double> CpAhconjHmconjVWp() const;
   double CphhconjHmconjVWp() const;
   double CpHmconjHmVP() const;
   double CpHmconjHmVZ() const;
   double CpHmconjHmconjVWpVWp() const;
   std::complex<double> CpHmconjHmVZVZ() const;
   std::complex<double> CpbarFuFdconjHmPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFdconjHmPL(int gI1, int gI2) const;
   std::complex<double> CpbarFvFeconjHmPR(int gI1, int gI2) const;
   double CpbarFvFeconjHmPL(int , int ) const;
   double CpAhAhhh() const;
   std::complex<double> CpbargWpgWpAh() const;
   std::complex<double> CpbargWpCgWpCAh() const;
   double CpAhAhAhAh() const;
   double CpAhAhhhhh() const;
   std::complex<double> CpAhhhVZ() const;
   std::complex<double> CpAhHmVWp() const;
   double CpAhAhconjVWpVWp() const;
   std::complex<double> CpAhAhVZVZ() const;
   std::complex<double> CpbarFdFdAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gI1, int gI2) const;
   double Cphhhhhh() const;
   double CphhVZVZ() const;
   double CphhconjVWpVWp() const;
   double CpbargWpgWphh() const;
   double CpbargWpCgWpChh() const;
   double CpbargZgZhh() const;
   double Cphhhhhhhh() const;
   double CphhHmVWp() const;
   double CphhhhconjVWpVWp() const;
   std::complex<double> CphhhhVZVZ() const;
   std::complex<double> CpbarFdFdhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gI1, int gI2) const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpHmVPVWp() const;
   double CpbargWpgWpVP() const;
   double CpbargWpCgWpCVP() const;
   std::complex<double> CpHmconjHmVPVP() const;
   double CpconjVWpVPVWp() const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   double CpconjVWpVPVPVWp3() const;
   double CpconjVWpVPVPVWp1() const;
   double CpconjVWpVPVPVWp2() const;
   double CpHmVWpVZ() const;
   double CpbargWpgWpVZ() const;
   double CpbargWpCgWpCVZ() const;
   double CpconjVWpVWpVZ() const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int , int ) const;
   double CpconjVWpVWpVZVZ1() const;
   double CpconjVWpVWpVZVZ2() const;
   double CpconjVWpVWpVZVZ3() const;
   double CpbargPgWpconjVWp() const;
   double CpbargWpCgPconjVWp() const;
   double CpbargWpCgZconjVWp() const;
   double CpbargZgWpconjVWp() const;
   std::complex<double> CpbarFdFuconjVWpPL(int gI1, int gI2) const;
   double CpbarFdFuconjVWpPR(int , int ) const;
   std::complex<double> CpbarFeFvconjVWpPL(int gI1, int gI2) const;
   double CpbarFeFvconjVWpPR(int , int ) const;
   double CpconjVWpconjVWpVWpVWp2() const;
   double CpconjVWpconjVWpVWpVWp1() const;
   double CpconjVWpconjVWpVWpVWp3() const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFuHmPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFuHmPR(int gO1, int gI2) const;
   double CpbarUFdFuconjVWpPR(int , int ) const;
   std::complex<double> CpbarUFdFuconjVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1) const;
   double CpbarUFuFdVWpPR(int , int ) const;
   std::complex<double> CpbarUFuFdVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFdconjHmPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFdconjHmPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFvHmPL(int gO2, int gI2) const;
   double CpbarUFeFvHmPR(int , int ) const;
   double CpbarUFeFvconjVWpPR(int , int ) const;
   double CpbarUFeFvconjVWpPL(int gO1, int gI2) const;
   double CpbarFvFeVWpPR(int , int ) const;
   std::complex<double> CpbarFvFeVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarFdFuHmPL(int gO2, int gI2) const;
   std::complex<double> CpbarFdFuHmPR(int gO1, int gI2) const;
   std::complex<double> CpbarFeFvHmPL(int gO2, int gI2) const;
   double CpbarFeFvHmPR(int , int ) const;
   double CpbarFuFdVWpPR(int , int ) const;
   std::complex<double> CpbarFuFdVWpPL(int gO1, int gI2) const;
   std::complex<double> CpFvFvUhhPL(int gI1, int gI2) const;
   std::complex<double> CpFvFvUhhPR(int gI1, int gI2) const;
   std::complex<double> CpFvFvUAhPL(int gI1, int gI2) const;
   std::complex<double> CpFvFvUAhPR(int gI1, int gI2) const;
   std::complex<double> self_energy_Hm_1loop(double p ) const;
   std::complex<double> self_energy_Ah_1loop(double p ) const;
   std::complex<double> self_energy_hh_1loop(double p ) const;
   std::complex<double> self_energy_VG_1loop(double p ) const;
   std::complex<double> self_energy_VP_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VWp_1loop(double p ) const;
   std::complex<double> self_energy_Fd_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL(double p) const;
   std::complex<double> self_energy_Fu_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL(double p) const;
   std::complex<double> self_energy_Fe_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL(double p) const;
   std::complex<double> self_energy_Fv_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_1(double p) const;
   std::complex<double> self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PR(double p) const;
   std::complex<double> self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PL(double p) const;
   std::complex<double> self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy(double p) const;
   std::complex<double> tadpole_hh_1loop() const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations() const;
   /// calculates the tadpoles divided by VEVs at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations_over_vevs() const;







   void calculate_MVG_pole();
   void calculate_MFv_pole();
   void calculate_Mhh_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_MVWp_pole();
   double calculate_MVWp_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double);
   double calculate_MVZ_DRbar(double);
   double calculate_MVWp_DRbar(double);

   double ThetaW() const;


private:
   int ewsb_loop_order{4};           ///< loop order for EWSB
   int pole_mass_loop_order{4};      ///< loop order for pole masses
   bool calculate_sm_pole_masses{false};  ///< switch to calculate the pole masses of the Standard Model particles
   bool calculate_bsm_pole_masses{true};  ///< switch to calculate the pole masses of the BSM particles
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< RG running precision
   double ewsb_iteration_precision{1.e-5};///< precision goal of EWSB solution
   cSMHdCKM_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{cSMHdCKM_info::model_name,
                     &cSMHdCKM_info::particle_names_getter,
                     &cSMHdCKM_info::parameter_names_getter}; ///< problems
   Loop_corrections loop_corrections{}; ///< used pole mass corrections
   std::shared_ptr<cSMHdCKM_ewsb_solver_interface> ewsb_solver{};
   Threshold_corrections threshold_corrections{}; ///< used threshold corrections

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_tree_level_custom();
   void copy_DRbar_masses_to_pole_masses();

   // Passarino-Veltman loop functions
   double A0(double) const noexcept;
   double B0(double, double, double) const noexcept;
   double B1(double, double, double) const noexcept;
   double B00(double, double, double) const noexcept;
   double B22(double, double, double) const noexcept;
   double H0(double, double, double) const noexcept;
   double F0(double, double, double) const noexcept;
   double G0(double, double, double) const noexcept;

   // DR-bar masses
   double MVG{};
   double MHm{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MAh{};
   double Mhh{};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   double MVWp{};
   double MVP{};
   double MVZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> UV{Eigen::Matrix<std::complex<double>,3,3>::Zero()};

   // phases

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const cSMHdCKM_mass_eigenstates&);

} // namespace flexiblesusy

#endif

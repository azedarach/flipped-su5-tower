#include "cSMHdCKMRHNEFT_matching.hpp"
#include "cSMHdCKM_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_info.hpp"

#include "config.h"
#include "global_thread_pool.hpp"
#include "linalg2.hpp"
#include "loop_corrections.hpp"
#include "single_scale_matching.hpp"
#include "wrappers.hpp"

#include <cmath>

namespace flexiblesusy {
namespace cSMHdCKMRHNEFT_matching {

#define MODELPARAMETER(p) model.get_##p()
#define EXTRAPARAMETER(parameter) model.get_##parameter()
#define EFTPARAMETER(p) eft.get_##p()
#define INPUTPARAMETER(p) model.get_input().p
#define PHASE(p) model.get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define Pole(p) model.get_physical().p
#define SCALE model.get_scale()

namespace {

/**
 * Returns EFT parameters with tree-level lambda.  The tree-level
 * lambda is calculated from the tree-level cSMHdCKMRHN parameters.
 *
 * @param eft EFT parameters
 * @param model_0l cSMHdCKMRHN tree-level parameters
 *
 * @return EFT tree-level parameters
 */
cSMHdCKM_mass_eigenstates calculate_EFT_tree_level(
   const cSMHdCKM_mass_eigenstates& eft,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   auto eft_0l = eft;
   eft_0l.get_problems().clear();
   eft_0l.set_Lambdax(Sqr(model_0l.get_Mhh()/eft.get_v()));
   eft_0l.set_pole_mass_loop_order(0);
   eft_0l.set_ewsb_loop_order(0);
   eft_0l.solve_ewsb_tree_level();
   eft_0l.calculate_DRbar_masses();

   return eft_0l;
}

/**
 * Calculates cSMHdCKMRHN tree-level parameters (g1, g2, g3, Yu, Yd,
 * Ye, v) by performing a tree-level matching from the EFT to the
 * cSMHdCKMRHN.
 *
 * @param model cSMHdCKMRHN parameters
 * @param eft EFT parameters
 *
 * @return cSMHdCKMRHN with tree-level parameters
 */
cSMHdCKMRHN_mass_eigenstates calculate_cSMHdCKMRHN_tree_level(
   const cSMHdCKMRHN_mass_eigenstates& model,
   const cSMHdCKM_mass_eigenstates& eft)
{
   auto model_0l = model;
   model_0l.get_problems().clear();
   model_0l.set_pole_mass_loop_order(0);
   model_0l.set_ewsb_loop_order(0);
   model_0l.solve_ewsb_tree_level();
   model_0l.calculate_DRbar_masses();

   match_low_to_high_scale_model_tree_level(model_0l, eft);

   return model_0l;
}

/**
 * Calculates squared Higgs pole mass in the EFT,
 * \f$(M_h^{\text{EFT}})^2\f$.
 *
 * @param eft_0l EFT parameters with tree-level lambda
 *
 * @return squared Higgs pole mass in the EFT
 */
double calculate_Mh2_pole(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double p = eft_0l.get_Mhh();
   const double self_energy = Re(eft_0l.self_energy_hh_1loop(p));
   const double tadpole = Re(eft_0l.tadpole_hh_1loop() / eft_0l.get_v());
   const double mh2_tree = Sqr(eft_0l.get_Mhh());
   const double Mh2_pole = mh2_tree - self_energy + tadpole;

   return Mh2_pole;
}

/**
 * Calculates tadpoles over vevs (at given fixed loop order)
 *
 * @param model model parameters
 * @param loop_order loop order
 *
 * @return tadpole at loop order
 */
Eigen::Matrix<double,1,1> calculate_tadpole_over_vevs(
   cSMHdCKMRHN_mass_eigenstates model, int loop_order)
{
   if (loop_order == 0) {
      model.set_ewsb_loop_order(loop_order);
      return model.tadpole_equations_over_vevs();
   }

   model.set_ewsb_loop_order(loop_order);
   const auto tadpole_lo = model.tadpole_equations_over_vevs(); // nL

   model.set_ewsb_loop_order(loop_order - 1);
   const auto tadpole_lom1 = model.tadpole_equations_over_vevs(); // (n-1)L

   return (tadpole_lo - tadpole_lom1).eval();
}

/**
 * Calculates tree-level Higgs mass matrix in the cSMHdCKMRHN without
 * including tadpoles implicitly.  This is achieved by solving the
 * EWSB equations at tree-level in order to avoid the inclusion of
 * loop tadpoles.
 *
 * @param model cSMHdCKMRHN parameters
 * @return tree-level Higgs mass matrix in the cSMHdCKMRHN w/o tadpoles
 */
auto calculate_mh2_tree_level(cSMHdCKMRHN_mass_eigenstates model) -> decltype(model.get_mass_matrix_hh())
{
   model.solve_ewsb_tree_level();
   return model.get_mass_matrix_hh();
}

/**
 * Calculates squared Higgs pole mass in the cSMHdCKMRHN,
 * \f$(M_h^{\text{cSMHdCKMRHN}})^2\f$.
 *
 * @param model_0l tree-level cSMHdCKMRHN parameters
 * @param model_1l 1-loop cSMHdCKMRHN parameters
 *
 * @return squared Higgs pole mass in the cSMHdCKMRHN
 */
double calculate_Mh2_pole(
   const cSMHdCKMRHN_mass_eigenstates& model_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_1l)
{
   // calculate tree-level mass matrix
   const auto mh2_tree = calculate_mh2_tree_level(model_1l);

   // calculate 1L self-energy using tree-level parameters
   const double p = model_0l.get_Mhh();
   const auto self_energy = Re(model_0l.self_energy_hh_1loop(p));

   // calculate 1L tadpoles using tree-level parameters
   const auto tadpole = calculate_tadpole_over_vevs(model_0l, 1);

   double Mh2_pole = 0.;

   Mh2_pole = mh2_tree - self_energy - tadpole(0);

   return Mh2_pole;
}

/**
 * Calculates \f$\lambda(Q)\f$ at the current loop level from the
 * lightest CP-even Higgs boson mass of the cSMHdCKMRHN by requiring
 * that the Higgs pole masses are equal in both models.
 *
 * @param eft EFT
 * @param model_1l cSMHdCKMRHN parameters
 */
void match_high_to_low_scale_model_1loop(
   cSMHdCKM_mass_eigenstates& eft,
   const cSMHdCKMRHN_mass_eigenstates& model_1l)
{
   // tree-level cSMHdCKMRHN parameters
   const auto model_0l = calculate_cSMHdCKMRHN_tree_level(model_1l, eft);
   const auto eft_0l = calculate_EFT_tree_level(eft, model_0l);

   const double mh2_eft = Sqr(eft_0l.get_Mhh());
   const double Mh2_eft = calculate_Mh2_pole(eft_0l);
   const double Mh2_bsm = calculate_Mh2_pole(model_0l, model_1l);

   eft.set_Lambdax((Mh2_bsm - Mh2_eft + mh2_eft)/Sqr(eft.get_v()));

   eft.get_problems().add(eft_0l.get_problems());
}

double calculate_delta_alpha_em(double alpha_em, const cSMHdCKMRHN_mass_eigenstates& model)
{
   const double currentScale = model.get_scale();
   double delta_alpha_em = 0.;

   
   delta_alpha_em += alpha_em/(2.*Pi)*(0);


   return delta_alpha_em;
}

double calculate_delta_alpha_s(double alpha_s, const cSMHdCKMRHN_mass_eigenstates& model)
{
   const double currentScale = model.get_scale();
   double delta_alpha_s = 0.;

   
   delta_alpha_s += alpha_s/(2.*Pi)*(0);


   return delta_alpha_s;
}

Eigen::Matrix<double,3,3> calculate_MFu_DRbar_tree_level(const cSMHdCKMRHN_mass_eigenstates& model)
{
   Eigen::Matrix<double,3,3> mf = ZEROMATRIX(3,3);

   mf(0,0) = model.get_MFu(0);
   mf(1,1) = model.get_MFu(1);
   mf(2,2) = model.get_MFu(2);


   return mf;
}

Eigen::Matrix<double,3,3> calculate_MFd_DRbar_tree_level(const cSMHdCKMRHN_mass_eigenstates& model)
{
   Eigen::Matrix<double,3,3> mf = ZEROMATRIX(3,3);

   mf(0,0) = model.get_MFd(0);
   mf(1,1) = model.get_MFd(1);
   mf(2,2) = model.get_MFd(2);


   return mf;
}

Eigen::Matrix<double,3,3> calculate_MFe_DRbar_tree_level(const cSMHdCKMRHN_mass_eigenstates& model)
{
   Eigen::Matrix<double,3,3> mf = ZEROMATRIX(3,3);

   mf(0,0) = model.get_MFe(0);
   mf(1,1) = model.get_MFe(1);
   mf(2,2) = model.get_MFe(2);


   return mf;
}

double calculate_MFu_pole_1loop(
   int i,
   const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double p = eft_0l.get_MFu(i);
   const auto self_energy_1  = Re(eft_0l.self_energy_Fu_1loop_1(p));
   const auto self_energy_PL = Re(eft_0l.self_energy_Fu_1loop_PL(p));
   const auto self_energy_PR = Re(eft_0l.self_energy_Fu_1loop_PR(p));
   const auto M_tree = eft_0l.get_mass_matrix_Fu();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> MFu_pole;
   fs_svd(M_loop, MFu_pole);

   return MFu_pole(i);
}

double calculate_MFd_pole_1loop(
   int i,
   const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double p = eft_0l.get_MFd(i);
   const auto self_energy_1  = Re(eft_0l.self_energy_Fd_1loop_1(p));
   const auto self_energy_PL = Re(eft_0l.self_energy_Fd_1loop_PL(p));
   const auto self_energy_PR = Re(eft_0l.self_energy_Fd_1loop_PR(p));
   const auto M_tree = eft_0l.get_mass_matrix_Fd();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> MFd_pole;
   fs_svd(M_loop, MFd_pole);

   return MFd_pole(i);
}

double calculate_MFe_pole_1loop(
   int i,
   const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double p = eft_0l.get_MFe(i);
   const auto self_energy_1  = Re(eft_0l.self_energy_Fe_1loop_1(p));
   const auto self_energy_PL = Re(eft_0l.self_energy_Fe_1loop_PL(p));
   const auto self_energy_PR = Re(eft_0l.self_energy_Fe_1loop_PR(p));
   const auto M_tree = eft_0l.get_mass_matrix_Fe();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> MFe_pole;
   fs_svd(M_loop, MFe_pole);

   return MFe_pole(i);
}

double calculate_MFu_pole_1loop(
   int i,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   double m_pole = 0.;

   const double p = model_0l.get_MFu(i);
   const auto self_energy_1  = Re(model_0l.self_energy_Fu_1loop_1(p));
   const auto self_energy_PL = Re(model_0l.self_energy_Fu_1loop_PL(p));
   const auto self_energy_PR = Re(model_0l.self_energy_Fu_1loop_PR(p));
   const auto M_tree = model_0l.get_mass_matrix_Fu();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> M_pole;
   fs_svd(M_loop, M_pole);

   m_pole = M_pole(i);

   return m_pole;
}

double calculate_MFd_pole_1loop(
   int i,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   double m_pole = 0.;

   const double p = model_0l.get_MFd(i);
   const auto self_energy_1  = Re(model_0l.self_energy_Fd_1loop_1(p));
   const auto self_energy_PL = Re(model_0l.self_energy_Fd_1loop_PL(p));
   const auto self_energy_PR = Re(model_0l.self_energy_Fd_1loop_PR(p));
   const auto M_tree = model_0l.get_mass_matrix_Fd();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> M_pole;
   fs_svd(M_loop, M_pole);

   m_pole = M_pole(i);

   return m_pole;
}

double calculate_MFe_pole_1loop(
   int i,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   double m_pole = 0.;

   const double p = model_0l.get_MFe(i);
   const auto self_energy_1  = Re(model_0l.self_energy_Fe_1loop_1(p));
   const auto self_energy_PL = Re(model_0l.self_energy_Fe_1loop_PL(p));
   const auto self_energy_PR = Re(model_0l.self_energy_Fe_1loop_PR(p));
   const auto M_tree = model_0l.get_mass_matrix_Fe();
   const auto M_loop = (M_tree - self_energy_PR * M_tree
                        - M_tree * self_energy_PL - self_energy_1).eval();

   Eigen::Array<double,3,1> M_pole;
   fs_svd(M_loop, M_pole);

   m_pole = M_pole(i);

   return m_pole;
}

Eigen::Matrix<double,3,3> calculate_MFu_pole_1loop(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFu_pole_1loop(0, eft_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFu_pole_1loop(1, eft_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFu_pole_1loop(2, eft_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFu_pole_1loop(0, eft_0l),
             calculate_MFu_pole_1loop(1, eft_0l),
             calculate_MFu_pole_1loop(2, eft_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFd_pole_1loop(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFd_pole_1loop(0, eft_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFd_pole_1loop(1, eft_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFd_pole_1loop(2, eft_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFd_pole_1loop(0, eft_0l),
             calculate_MFd_pole_1loop(1, eft_0l),
             calculate_MFd_pole_1loop(2, eft_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFe_pole_1loop(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFe_pole_1loop(0, eft_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFe_pole_1loop(1, eft_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&eft_0l]{ return calculate_MFe_pole_1loop(2, eft_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFe_pole_1loop(0, eft_0l),
             calculate_MFe_pole_1loop(1, eft_0l),
             calculate_MFe_pole_1loop(2, eft_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFu_pole_1loop(
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFu_pole_1loop(0, model_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFu_pole_1loop(1, model_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFu_pole_1loop(2, model_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFu_pole_1loop(0, model_0l),
             calculate_MFu_pole_1loop(1, model_0l),
             calculate_MFu_pole_1loop(2, model_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFd_pole_1loop(
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFd_pole_1loop(0, model_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFd_pole_1loop(1, model_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFd_pole_1loop(2, model_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFd_pole_1loop(0, model_0l),
             calculate_MFd_pole_1loop(1, model_0l),
             calculate_MFd_pole_1loop(2, model_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFe_pole_1loop(
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,1> M_pole;

#ifdef ENABLE_THREADS
   auto M_0 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFe_pole_1loop(0, model_0l); });
   auto M_1 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFe_pole_1loop(1, model_0l); });
   auto M_2 = global_thread_pool().run_packaged_task([&model_0l]{ return calculate_MFe_pole_1loop(2, model_0l); });
   M_pole << M_0.get(), M_1.get(), M_2.get();
#else
   M_pole << calculate_MFe_pole_1loop(0, model_0l),
             calculate_MFe_pole_1loop(1, model_0l),
             calculate_MFe_pole_1loop(2, model_0l);
#endif

   return M_pole.asDiagonal();
}

Eigen::Matrix<double,3,3> calculate_MFu_DRbar_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,3> mf_eft = ZEROMATRIX(3,3);

   const auto Mf_eft  = calculate_MFu_pole_1loop(eft_0l);
   const auto Mf_bsm = calculate_MFu_pole_1loop(model_0l);
   const auto mf_bsm = calculate_MFu_DRbar_tree_level(model_0l);

   mf_eft = mf_bsm - Mf_bsm + Mf_eft;

   return Abs(mf_eft);
}

Eigen::Matrix<double,3,3> calculate_MFd_DRbar_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,3> mf_eft = ZEROMATRIX(3,3);

   const auto Mf_eft  = calculate_MFd_pole_1loop(eft_0l);
   const auto Mf_bsm = calculate_MFd_pole_1loop(model_0l);
   const auto mf_bsm = calculate_MFd_DRbar_tree_level(model_0l);

   mf_eft = mf_bsm - Mf_bsm + Mf_eft;

   return Abs(mf_eft);
}

Eigen::Matrix<double,3,3> calculate_MFe_DRbar_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   Eigen::Matrix<double,3,3> mf_eft = ZEROMATRIX(3,3);

   const auto Mf_eft  = calculate_MFe_pole_1loop(eft_0l);
   const auto Mf_bsm = calculate_MFe_pole_1loop(model_0l);
   const auto mf_bsm = calculate_MFe_DRbar_tree_level(model_0l);

   mf_eft = mf_bsm - Mf_bsm + Mf_eft;

   return Abs(mf_eft);
}

double calculate_MW_pole_1loop(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double mw = eft_0l.get_MVWp();
   const double p = eft_0l.get_MVWp();
   const double self_energy = Re(eft_0l.self_energy_VWp_1loop(p));
   const double M_loop = Sqr(mw) - self_energy;

   return SignedAbsSqrt(M_loop);
}

double calculate_MW_pole_1loop(const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   const double mw = model_0l.get_MVWp();
   const double p = model_0l.get_MVWp();
   const double self_energy = Re(model_0l.self_energy_VWp_1loop(p));
   const double M_loop = Sqr(mw) - self_energy;

   return SignedAbsSqrt(M_loop);
}

double calculate_MZ_pole_1loop(const cSMHdCKM_mass_eigenstates& eft_0l)
{
   const double mz = eft_0l.get_MVZ();
   const double p = eft_0l.get_MVZ();
   const double self_energy = Re(eft_0l.self_energy_VZ_1loop(p));
   const double M_loop = Sqr(mz) - self_energy;

   return SignedAbsSqrt(M_loop);
}

double calculate_MZ_pole_1loop(const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   const double mz = model_0l.get_MVZ();
   const double p = model_0l.get_MVZ();
   const double self_energy = Re(model_0l.self_energy_VZ_1loop(p));
   const double M_loop = Sqr(mz) - self_energy;

   return SignedAbsSqrt(M_loop);
}

double calculate_MW_DRbar_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   const double MW_eft = calculate_MW_pole_1loop(eft_0l);
   const double MW_bsm = calculate_MW_pole_1loop(model_0l);
   const double mw2 = Sqr(MW_eft) - Sqr(MW_bsm) + Sqr(model_0l.get_MVWp());

   return AbsSqrt(mw2);
}

double calculate_MZ_DRbar_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l)
{
   const double MZ_eft = calculate_MZ_pole_1loop(eft_0l);
   const double MZ_bsm = calculate_MZ_pole_1loop(model_0l);
   const double mz2 = Sqr(MZ_eft) - Sqr(MZ_bsm) + Sqr(model_0l.get_MVZ());

   return AbsSqrt(mz2);
}

/**
 * Calculates cSMHdCKMRHN parameters at 1-loop level by performing a
 * 1-loop matching.
 *
 * @param eft_0l EFT parameters (lambda is at tree-level)
 * @param eft_1l EFT parameters (level is at 1-loop level)
 * @param model_0l cSMHdCKMRHN parameters at tree-level
 * @param model_1l cSMHdCKMRHN parameters at 1-loop level
 *
 * @return cSMHdCKMRHN 1-loop parameters
 */
cSMHdCKMRHN_mass_eigenstates calculate_cSMHdCKMRHN_1loop(
   const cSMHdCKM_mass_eigenstates& eft_0l,
   const cSMHdCKM_mass_eigenstates& eft_1l,
   const cSMHdCKMRHN_mass_eigenstates& model_0l,
   const cSMHdCKMRHN_mass_eigenstates& model_1l)
{
   auto model = model_1l;
   auto eft = eft_1l;

   model.calculate_DRbar_masses();
   model.solve_ewsb();

   eft.calculate_DRbar_masses();
   eft.solve_ewsb();

   const double alpha_em = Sqr(eft_0l.get_g1() * eft_0l.get_g2() * cSMHdCKM_info::normalization_g1 * cSMHdCKM_info::normalization_g2)
            /(4. * Pi * (Sqr(eft_0l.get_g1()*cSMHdCKM_info::normalization_g1) + Sqr(eft_0l.get_g2()*cSMHdCKM_info::normalization_g2)));
   const double alpha_s  = Sqr(eft_0l.get_g3() * cSMHdCKM_info::normalization_g3)/(4. * Pi);
   const double delta_alpha_em = calculate_delta_alpha_em(alpha_em, model_0l);
   const double delta_alpha_s = calculate_delta_alpha_s(alpha_s, model_0l);

   // running cSMHdCKMRHN W, Z masses (via 1L matching)
   const double mW2_1L = Sqr(calculate_MW_DRbar_1loop(eft_0l, model_0l));
   const double mZ2_1L = Sqr(calculate_MZ_DRbar_1loop(eft_0l, model_0l));

   // running cSMHdCKMRHN quark and lepton masses (via 1L matching)
   const Eigen::Matrix<double, 3, 3> upQuarksDRbar    = calculate_MFu_DRbar_1loop(eft_0l, model_0l);
   const Eigen::Matrix<double, 3, 3> downQuarksDRbar  = calculate_MFd_DRbar_1loop(eft_0l, model_0l);
   const Eigen::Matrix<double, 3, 3> downLeptonsDRbar = calculate_MFe_DRbar_1loop(eft_0l, model_0l);

   // running cSMHdCKMRHN gauge couplings (via 1L matching)
   const double g1_1L = AbsSqrt(4. * Pi * alpha_em * (1. + delta_alpha_em) * mZ2_1L / mW2_1L) / cSMHdCKMRHN_info::normalization_g1;
   const double g2_1L = AbsSqrt(4. * Pi * alpha_em * (1. + delta_alpha_em) / (1. - mW2_1L/mZ2_1L)) / cSMHdCKMRHN_info::normalization_g2;
   const double g3_1L = AbsSqrt(4. * Pi * alpha_s * (1. + delta_alpha_s)) / cSMHdCKMRHN_info::normalization_g3;

   model.set_g1(g1_1L);
   model.set_g2(g2_1L);
   model.set_g3(g3_1L);

   {
      auto MODEL = &model;
      const double VEV = 2. * AbsSqrt(mZ2_1L/(Sqr(g1_1L*cSMHdCKMRHN_info::normalization_g1) + Sqr(g2_1L*cSMHdCKMRHN_info::normalization_g2)));

      
      MODEL->set_v(Re(VEV));

   }

   const auto v = MODELPARAMETER(v);
   model.set_Yu(((1.4142135623730951*upQuarksDRbar)/v).transpose());
   model.set_Yd(((1.4142135623730951*downQuarksDRbar)/v).transpose());
   model.set_Ye(((1.4142135623730951*downLeptonsDRbar)/v).transpose());


   return model;
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the cSMHdCKMRHN at the necessary loop level from the known Standard
 * Model couplings and the SM vev.
 *
 * @param model cSMHdCKMRHN parameters (to be set)
 * @param eft EFT parameters
 */
void match_low_to_high_scale_model_1loop(
   cSMHdCKMRHN_mass_eigenstates& model,
   const cSMHdCKM_mass_eigenstates& eft)
{
   // tree-level cSMHdCKMRHN parameters
   const auto model_0l = calculate_cSMHdCKMRHN_tree_level(model, eft);
   const auto eft_0l = calculate_EFT_tree_level(eft, model_0l);

   // 1-loop parameters
   const auto model_1l = calculate_cSMHdCKMRHN_1loop(eft_0l, eft, model_0l, model);

   model.set_g1(model_0l.get_g1());
   model.set_g2(model_0l.get_g2());
   model.set_g3(model_0l.get_g3());
   model.set_v(model_1l.get_v());
   model.set_Yd(model_0l.get_Yd());
   model.set_Ye(model_0l.get_Ye());
   model.set_Yu(model_0l.get_Yu());


   model.get_problems().add(model_0l.get_problems());
   model.get_problems().add(model_1l.get_problems());
}

} // anonymous namespace

/**
 * Calculates \f$\lambda(Q)\f$ at the tree level from the lightest
 * CP-even Higgs boson mass of the cSMHdCKMRHN.
 *
 * @param eft EFT parameters
 * @param model cSMHdCKMRHN parameters
 */
void match_high_to_low_scale_model_tree_level(
   cSMHdCKM_mass_eigenstates& eft, const cSMHdCKMRHN_mass_eigenstates& model)
{
   auto model_tmp = model;
   model_tmp.calculate_DRbar_masses();
   eft.set_Lambdax(Sqr(model_tmp.get_Mhh()/eft.get_v()));
   eft.calculate_DRbar_masses();
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the cSMHdCKMRHN at the tree level from the known Standard Model
 * couplings and the SM vev.
 */
void match_low_to_high_scale_model_tree_level(
   cSMHdCKMRHN_mass_eigenstates& model, const cSMHdCKM_mass_eigenstates& eft_input)
{
   auto eft = eft_input;
   eft.calculate_DRbar_masses();

   model.set_g1(eft.get_g1()*cSMHdCKM_info::normalization_g1/cSMHdCKMRHN_info::normalization_g1);
   model.set_g2(eft.get_g2()*cSMHdCKM_info::normalization_g2/cSMHdCKMRHN_info::normalization_g2);
   model.set_g3(eft.get_g3()*cSMHdCKM_info::normalization_g3/cSMHdCKMRHN_info::normalization_g3);

   {
      auto MODEL = &model;
      const double VEV = eft.get_v();

      
      MODEL->set_v(Re(VEV));

   }

   Eigen::Matrix<double, 3, 3> upQuarksDRbar    = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downQuarksDRbar  = ZEROMATRIX(3,3);
   Eigen::Matrix<double, 3, 3> downLeptonsDRbar = ZEROMATRIX(3,3);

   upQuarksDRbar.diagonal()    = eft.get_MFu();
   downQuarksDRbar.diagonal()  = eft.get_MFd();
   downLeptonsDRbar.diagonal() = eft.get_MFe();

   const auto v = MODELPARAMETER(v);
   model.set_Yu(((1.4142135623730951*upQuarksDRbar)/v).transpose());
   model.set_Yd(((1.4142135623730951*downQuarksDRbar)/v).transpose());
   model.set_Ye(((1.4142135623730951*downLeptonsDRbar)/v).transpose());


   model.calculate_DRbar_masses();
}

/**
 * Calculates the gauge and Yukawa couplings and the SM-like VEV in
 * the cSMHdCKMRHN at the 1-loop level from the known Standard Model
 * couplings and the SM vev.
 *
 * @note The \a loop_order argument is interpreted to be the loop
 * order of the "downwards matching".  The functions responsible for
 * the upwards matching will use the downwards matching loop order to
 * determine the cSMHdCKMRHN at the right loop order to avoid the
 * appearence of large logarithms.
 *
 * @param model cSMHdCKMRHN parameters
 * @param eft EFT
 * @param loop_order downwards matching loop order
 */
void match_low_to_high_scale_model(
   cSMHdCKMRHN_mass_eigenstates& model, const cSMHdCKM_mass_eigenstates& eft, int loop_order)
{
   if (loop_order == 0) {
      match_low_to_high_scale_model_tree_level(model, eft);
      return;
   }

   match_low_to_high_scale_model_1loop(model, eft);
}

/**
 * Calculates \f$\lambda(Q)\f$ at the 1-loop level from the lightest
 * CP-even Higgs boson mass of the cSMHdCKMRHN by requiring that the
 * 1-loop Higgs pole masses are equal in both models.
 *
 * @param eft EFT
 * @param model cSMHdCKMRHN parameters
 * @param loop_order downwards matching loop order
 */
void match_high_to_low_scale_model(
   cSMHdCKM_mass_eigenstates& eft, const cSMHdCKMRHN_mass_eigenstates& model, int loop_order)
{
   if (loop_order == 0) {
      match_high_to_low_scale_model_tree_level(eft, model);
      return;
   }

   match_high_to_low_scale_model_1loop(eft, model);
}

} // namespace cSMHdCKMRHNEFT_matching
} // namespace flexiblesusy

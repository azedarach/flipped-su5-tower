#ifndef cSMHdCKMRHNEFT_TWO_SCALE_INITIAL_GUESSER_H
#define cSMHdCKMRHNEFT_TWO_SCALE_INITIAL_GUESSER_H

#include "cSMHdCKMRHN_initial_guesser.hpp"
#include "cSMHdCKMRHN_two_scale_susy_scale_constraint.hpp"
#include "cSMHdCKMRHN_two_scale_high_scale_constraint.hpp"
#include "cSMHdCKM_two_scale_low_scale_constraint.hpp"
#include "initial_guesser.hpp"
#include "lowe.h"

#include <sstream>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template <class T>
class cSMHdCKMRHN;

template <class T>
class cSMHdCKM;

template <class T>
class cSMHdCKMRHNEFT_initial_guesser;

/**
 * @class cSMHdCKMRHNEFT_initial_guesser<Two_scale>
 * @brief initial guesser for the cSMHdCKMRHNEFT tower
 */

template<>
class cSMHdCKMRHNEFT_initial_guesser<Two_scale> : public Initial_guesser {
public:
   cSMHdCKMRHNEFT_initial_guesser(cSMHdCKMRHN<Two_scale>*,
                                  cSMHdCKM<Two_scale>*,
                                  const softsusy::QedQcd&,
                                  const cSMHdCKM_low_scale_constraint<Two_scale>&,
                                  const cSMHdCKMRHN_susy_scale_constraint<Two_scale>&,
                                  const cSMHdCKMRHN_high_scale_constraint<Two_scale>&);
   virtual ~cSMHdCKMRHNEFT_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   cSMHdCKMRHN<Two_scale>* model{nullptr}; ///< pointer to model class
   cSMHdCKM<Two_scale>* eft{nullptr}; ///< pointer to effective low energy model
   softsusy::QedQcd qedqcd{}; ///< Standard Model low-energy data
   double mu_guess{0.}; ///< guessed DR-bar mass of up-quark
   double mc_guess{0.}; ///< guessed DR-bar mass of charm-quark
   double mt_guess{0.}; ///< guessed DR-bar mass of top-quark
   double md_guess{0.}; ///< guessed DR-bar mass of down-quark
   double ms_guess{0.}; ///< guessed DR-bar mass of strange-quark
   double mb_guess{0.}; ///< guessed DR-bar mass of bottom-quark
   double me_guess{0.}; ///< guessed DR-bar mass of electron
   double mm_guess{0.}; ///< guessed DR-bar mass of muon
   double mtau_guess{0.}; ///< guessed DR-bar mass of tau
   double running_precision{1.0e-3}; ///< Runge-Kutta RG running precision
   cSMHdCKM_low_scale_constraint<Two_scale> low_constraint{};
   cSMHdCKMRHN_susy_scale_constraint<Two_scale> susy_constraint{};
   cSMHdCKMRHN_high_scale_constraint<Two_scale> high_constraint{};
   Eigen::Matrix<std::complex<double>,3,3> upQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downQuarksDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> downLeptonsDRbar{Eigen::Matrix<std::complex<double>,3,3>::Zero()};

   void guess_eft_parameters();
   void guess_model_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_running_SM_masses();
};

} // namespace flexiblesusy

#endif

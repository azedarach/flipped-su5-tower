#include "cSMHdCKMRHNEFT_two_scale_initial_guesser.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"
#include "cSMHdCKMRHNEFT_two_scale_matching.hpp"
#include "cSMHdCKM_two_scale_model.hpp"

#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
// @todo remove
#include <iostream>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define SMPARAMETER(p) eft->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

cSMHdCKMRHNEFT_initial_guesser<Two_scale>::cSMHdCKMRHNEFT_initial_guesser(
   cSMHdCKMRHN<Two_scale>* model_,
   cSMHdCKM<Two_scale>* eft_,
   const softsusy::QedQcd& qedqcd_,
   const cSMHdCKMRHNEFT_low_scale_constraint<Two_scale>& low_constraint_,
   const cSMHdCKMRHNEFT_susy_scale_constraint<Two_scale>& susy_constraint_,
   const cSMHdCKMRHNEFT_high_scale_constraint<Two_scale>& high_constraint_
)
   : Initial_guesser()
   , model(model_)
   , eft(eft_)
   , qedqcd(qedqcd_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
{
   if (!model)
      throw SetupError("cSMHdCKMRHNEFT_initial_guesser: Error: pointer to model"
                       " cSMHdCKMRHN<Two_scale> must not be zero");
}

cSMHdCKMRHNEFT_initial_guesser<Two_scale>::~cSMHdCKMRHNEFT_initial_guesser()
{
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::guess()
{
   guess_eft_parameters();
   guess_model_parameters();
}

/**
 * Guesses the SUSY parameters (gauge, Yukawa couplings) at
 * \f$m_\text{top}^\text{pole}\f$ from the Standard Model gauge
 * couplings and fermion masses.  Threshold corrections are ignored.
 * The user-defined initial guess at the low-scale
 * (InitialGuessAtLowScale) is applied here:
 *
 * \code{.cpp}
   

 * \endcode
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::guess_eft_parameters()
{
   softsusy::QedQcd leAtMt(qedqcd);
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(softsusy::mUp);
   mc_guess = leAtMt.displayMass(softsusy::mCharm);
   mt_guess = model->get_thresholds() > 0 && model->get_threshold_corrections().mt > 0 ?
      leAtMt.displayMass(softsusy::mTop) - 30.0 :
      leAtMt.displayPoleMt();
   md_guess = leAtMt.displayMass(softsusy::mDown);
   ms_guess = leAtMt.displayMass(softsusy::mStrange);
   mb_guess = leAtMt.displayMass(softsusy::mBottom);
   me_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mElectron) :
      leAtMt.displayPoleMel();
   mm_guess = model->get_thresholds() > 0 ?
      leAtMt.displayMass(softsusy::mMuon) :
      leAtMt.displayPoleMmuon();
   mtau_guess = leAtMt.displayMass(softsusy::mTau);
   mv1_guess = leAtMt.displayNeutrinoPoleMass(1);
   mv2_guess = leAtMt.displayNeutrinoPoleMass(2);
   mv3_guess = leAtMt.displayNeutrinoPoleMass(3);

   calculate_running_SM_masses();

   // guess gauge couplings at mt
   const auto alpha_sm(leAtMt.guess_alpha_SM5(mtpole));

   eft->set_g1(sqrt(4.0 * Pi * alpha_sm(0)));
   eft->set_g2(sqrt(4.0 * Pi * alpha_sm(1)));
   eft->set_g3(sqrt(4.0 * Pi * alpha_sm(2)));
   eft->set_scale(mtpole);

   eft->set_v(Electroweak_constants::vev);
   eft->set_Yu(ZEROMATRIX(3,3));
   eft->set_Yd(ZEROMATRIX(3,3));
   eft->set_Ye(ZEROMATRIX(3,3));
   eft->set_Kappa(ZEROMATRIX(3,3));

   eft->set_Yu(0, 0, Sqrt(2.)* mu_guess/ eft->get_v());
   eft->set_Yu(1, 1, Sqrt(2.)* mc_guess/ eft->get_v());
   eft->set_Yu(2, 2, Sqrt(2.)* mt_guess/ eft->get_v());

   eft->set_Yd(0, 0, Sqrt(2.)* md_guess/ eft->get_v());
   eft->set_Yd(1, 1, Sqrt(2.)* ms_guess/ eft->get_v());
   eft->set_Yd(2, 2, Sqrt(2.)* mb_guess/ eft->get_v());

   eft->set_Ye(0, 0, Sqrt(2.)* me_guess/ eft->get_v());
   eft->set_Ye(1, 1, Sqrt(2.)* mm_guess/ eft->get_v());
   eft->set_Ye(2, 2, Sqrt(2.)* mtau_guess/ eft->get_v());

   eft->set_Kappa(0, 0, 4. * mv1_guess / Sqr(eft->get_v()));
   eft->set_Kappa(1, 1, 4. * mv2_guess / Sqr(eft->get_v()));
   eft->set_Kappa(2, 2, 4. * mv3_guess / Sqr(eft->get_v()));

   eft->set_Lambdax(0.12604);
   eft->solve_ewsb_tree_level();

   std::cout << "end of guess eft parameters:\n";
   eft->print(std::cout);
   std::cout << "check tree-level masses:\n";
   eft->calculate_MFu();
   eft->calculate_MFd();
   eft->calculate_MFe();
   eft->calculate_MFv();

   std::cout << "eft MFu = " << eft->get_MFu().transpose() << '\n';
   std::cout << "expected = " << mu_guess << ", " << mc_guess << ", " << mt_guess << '\n';
   std::cout << "eft MFd = " << eft->get_MFd().transpose() << '\n';
   std::cout << "expected = " << md_guess << ", " << ms_guess << ", " << mb_guess << '\n';
   std::cout << "eft MFe = " << eft->get_MFe().transpose() << '\n';
   std::cout << "expected = " << me_guess << ", " << mm_guess << ", " << mtau_guess << '\n';
   std::cout << "eft MFv = " << eft->get_MFv().transpose() << '\n';
   std::cout << "expected = " << mv1_guess << ", " << mv2_guess << ", " << mv3_guess << '\n';
}

void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_running_SM_masses();
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::calculate_running_SM_masses()
{
   upQuarksDRbar.setZero();
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;

   downQuarksDRbar.setZero();
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;

   downLeptonsDRbar.setZero();
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;

   neutrinosDRbar.setZero();
   neutrinosDRbar(0,0) = mv1_guess;
   neutrinosDRbar(1,1) = mv2_guess;
   neutrinosDRbar(2,2) = mv3_guess;
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::calculate_Yu_DRbar()
{

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::calculate_Yd_DRbar()
{

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::calculate_Ye_DRbar()
{

}

/**
 * Guesses the full model parameters.  At first it runs the effective SM to the
 * guess of the high-scale (SUSYScaleFirstGuess), and matches gauge and Yukawas
 * to the full theory. Then, initial guess and constrints are imposed at this scale:
 * \code{.cpp}
   
   MODEL->set_v(Re(LowEnergyConstant(vev)));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();

 * \endcode
 * After that, the model is run to the high scale, initial guess and constraints
 * are imposed (HighScaleInput):
 *
 * \code{.cpp}
   

 * \endcode
 *
 * Afterwards, it runs to the SUSY-scale guess (SUSYScaleFirstGuess) and
 * solves the EWSB conditions at the tree-level.  Finally the DR-bar
 * mass spectrum is calculated.
 */
void cSMHdCKMRHNEFT_initial_guesser<Two_scale>::guess_model_parameters()
{
   const double susy_scale_guess = susy_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_initial_scale_guess();
   const auto scale_getter = [this] () { return susy_constraint.get_scale(); };

   model->set_scale(susy_scale_guess);

   // apply susy-scale first guess
   {
   
   MODEL->set_v(Re(LowEnergyConstant(vev)));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_Yv(eft->get_Yu().transpose());
   }

   std::cout << "after apply susy-scale first guess:\n";
   model->print(std::cout);

   eft->run_to(susy_scale_guess, running_precision);
   eft->calculate_DRbar_masses();

   std::cout << "after running EFT to susy scale guess:\n";
   eft->print(std::cout);

   //get gauge and Yukawa couplings from effective theory
   std::cout << "performing tree-level matching\n";
   cSMHdCKMRHNEFT_matching_up<Two_scale> matching_up;
   matching_up.set_models(eft, model);
   matching_up.set_scale(scale_getter);
   matching_up.match_tree_level();

   std::cout << "after tree-level matching:\n";
   eft->print(std::cout);
   model->print(std::cout);

   model->run_to(susy_scale_guess, running_precision);

   // apply susy-scale constraint
   susy_constraint.set_model(model);
   susy_constraint.apply();

   std::cout << "after applying susy scale constraint:\n";
   model->print(std::cout);

   // run to high scale
   model->run_to(high_scale_guess, running_precision);

   // apply user-defined initial guess at the high scale
   {
   

   }

   // apply high-scale constraint
   high_constraint.set_model(model);
   high_constraint.apply();

   std::cout << "after applying high scale constraint:\n";
   model->print(std::cout);

   model->run_to(susy_scale_guess, running_precision);

   // apply EWSB constraint
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();

   std::cout << "end of guess model parameters:\n";
   model->print(std::cout);
}

} // namespace flexiblesusy

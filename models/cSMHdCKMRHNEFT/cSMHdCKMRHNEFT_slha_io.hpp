#ifndef cSMHdCKMRHNEFT_SLHA_IO_H
#define cSMHdCKMRHNEFT_SLHA_IO_H

#include "cSMHdCKMRHN_mass_eigenstates.hpp"
#include "cSMHdCKMRHN_model_slha.hpp"
#include "cSMHdCKMRHN_info.hpp"
#include "cSMHdCKMRHN_observables.hpp"
#include "cSMHdCKMRHN_physical.hpp"

#include "cSMHdCKM_model_slha.hpp"

#include "problems.hpp"
#include "spectrum_generator_problems.hpp"
#include "standard_model_two_scale_model.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "ew_input.hpp"
#include "lowe.h"

#include <Eigen/Core>
#include <string>
#include <tuple>
#include <utility>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODEL model
#define MODELPARAMETER(p) model.get_##p()
#define EXTRAPARAMETER(p) model.get_##p()
#define OBSERVABLES observables
#define LowEnergyConstant(p) Electroweak_constants::p
#define SCALES(p) scales.p

namespace flexiblesusy {

struct cSMHdCKMRHNEFT_input_parameters;
class Spectrum_generator_settings;

template <class T>
class cSMHdCKM;

template <class T>
class cSMHdCKMRHN;

struct cSMHdCKMRHNEFT_scales {
   double HighScale{0.}, SUSYScale{0.}, LowScale{0.};
   double pole_mass_scale{0.};
};

class cSMHdCKMRHNEFT_slha_io {
public:
   cSMHdCKMRHNEFT_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(cSMHdCKMRHNEFT_input_parameters&) const;
   void fill(cSMHdCKMRHN_mass_eigenstates&) const;
   template <class Model> void fill(cSMHdCKMRHN_slha<Model>&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   template <class Model> void set_extra(const cSMHdCKMRHN_slha<Model>&, const cSMHdCKMRHNEFT_scales&, const cSMHdCKMRHN_observables&);
   void set_input(const cSMHdCKMRHNEFT_input_parameters&);
   void set_modsel(const SLHA_io::Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class... Ts> void set_spectrum(const std::tuple<Ts...>&);
   template <class Model> void set_spectrum(const cSMHdCKM_slha<Model>&);
   template <class T> void set_spectrum(const cSMHdCKM<T>&);
   template <class Model> void set_spectrum(const cSMHdCKMRHN_slha<Model>&);
   template <class T> void set_spectrum(const cSMHdCKMRHN<T>&);
   void set_spectrum(const standard_model::Standard_model&);
   void set_spinfo(const Spectrum_generator_problems&);
   void set_spinfo(const Problems&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&) const;
   void write_to_file(const std::string& file_name) const { slha_io.write_to_file(file_name); }
   void write_to_stream(std::ostream& ostr = std::cout) const { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(cSMHdCKMRHNEFT_input_parameters&, int, double);
   static void fill_extpar_tuple(cSMHdCKMRHNEFT_input_parameters&, int, double);
   static void fill_imminpar_tuple(cSMHdCKMRHNEFT_input_parameters&, int, double);
   static void fill_imextpar_tuple(cSMHdCKMRHNEFT_input_parameters&, int, double);

   template <class EFT, class Model>
   static void fill_slhaea(SLHAea::Coll&, const cSMHdCKM_slha<EFT>&, const cSMHdCKMRHN_slha<Model>&,
                           const softsusy::QedQcd&, const cSMHdCKMRHNEFT_input_parameters&,
                           const cSMHdCKMRHNEFT_scales&, const cSMHdCKMRHN_observables&);

   template <class EFT, class Model>
   static SLHAea::Coll fill_slhaea(const cSMHdCKM_slha<EFT>&, const cSMHdCKMRHN_slha<Model>&,
                                   const softsusy::QedQcd&, const cSMHdCKMRHNEFT_input_parameters&,
                                   const cSMHdCKMRHNEFT_scales&, const cSMHdCKMRHN_observables&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;

   void set_extpar(const cSMHdCKMRHNEFT_input_parameters&);
   void set_imminpar(const cSMHdCKMRHNEFT_input_parameters&);
   void set_imextpar(const cSMHdCKMRHNEFT_input_parameters&);
   void set_minpar(const cSMHdCKMRHNEFT_input_parameters&);
   void set_mass(const cSMHdCKM_physical&, bool);
   void set_mass(const cSMHdCKMRHN_physical&, bool);
   void set_mass(const standard_model::Standard_model_physical&);
   void set_mixing_matrices(const cSMHdCKM_physical&, bool);
   void set_mixing_matrices(const cSMHdCKMRHN_physical&, bool);
   void set_mixing_matrices(const standard_model::Standard_model_physical&);
   template <class Model> void set_model_parameters(const cSMHdCKM_slha<Model>&);
   template <class Model> void set_model_parameters(const cSMHdCKMRHN_slha<Model>&);
   void set_model_parameters(const standard_model::Standard_model&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(cSMHdCKMRHN_mass_eigenstates&) const;
   void fill_physical(cSMHdCKMRHN_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::fill(cSMHdCKMRHN_slha<Model>& model) const
{
   fill(static_cast<cSMHdCKMRHN_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class EFT, class Model>
void cSMHdCKMRHNEFT_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const cSMHdCKM_slha<EFT>& eft,
   const cSMHdCKMRHN_slha<Model>& model,
   const softsusy::QedQcd& qedqcd, const cSMHdCKMRHNEFT_input_parameters& input,
   const cSMHdCKMRHNEFT_scales& scales, const cSMHdCKMRHN_observables& observables)
{
   cSMHdCKMRHNEFT_slha_io slha_io;
   const auto& problems = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_input(input);
   if (!error) {
      slha_io.set_spectrum(eft);
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales, observables);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class EFT, class Model>
SLHAea::Coll cSMHdCKMRHNEFT_slha_io::fill_slhaea(
   const cSMHdCKM_slha<EFT>& eft, const cSMHdCKMRHN_slha<Model>& model,
   const softsusy::QedQcd& qedqcd, const cSMHdCKMRHNEFT_input_parameters& input,
   const cSMHdCKMRHNEFT_scales& scales, const cSMHdCKMRHN_observables& observables)
{
   SLHAea::Coll slhaea;
   cSMHdCKMRHNEFT_slha_io::fill_slhaea(slhaea, eft, model, qedqcd, input, scales, observables);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::set_model_parameters(const cSMHdCKM_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block EFTgauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("EFTYu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block_imag("IMEFTYu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("EFTYd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block_imag("IMEFTYd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("EFTYe", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block_imag("IMEFTYe", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block EFTSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Lambdax)), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("EFTKappa", ToMatrix(MODELPARAMETER(Kappa)), "Kappa", model.get_scale());
   slha_io.set_block_imag("IMEFTKappa", ToMatrix(MODELPARAMETER(Kappa)), "Kappa", model.get_scale());
   {
      std::ostringstream block;
      block << "Block EFTHMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (MODELPARAMETER(v)), "v")
      ;
      slha_io.set_block(block);
   }


}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::set_model_parameters(const cSMHdCKMRHN_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block_imag("IMYu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block_imag("IMYd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block_imag("IMYe", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block("Yv", MODELPARAMETER(Yv), "Yv", model.get_scale());
   slha_io.set_block_imag("IMYv", MODELPARAMETER(Yv), "Yv", model.get_scale());
   slha_io.set_block("Mv", MODELPARAMETER(Mv), "Mv", model.get_scale());
   slha_io.set_block_imag("IMMv", MODELPARAMETER(Mv), "Mv", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Lambdax)), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (MODELPARAMETER(v)), "v")
      ;
      slha_io.set_block(block);
   }


}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 * @param scales struct of boundary condition scales
 * @param observables struct of observables
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::set_extra(
   const cSMHdCKMRHN_slha<Model>& model, const cSMHdCKMRHNEFT_scales& scales,
   const cSMHdCKMRHN_observables& observables)
{
   const cSMHdCKMRHN_physical physical(model.get_physical_slha());

   {
      std::ostringstream block;
      block << "Block FlexibleSUSYLowEnergy Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(21, (OBSERVABLES.a_muon), "Delta(g-2)_muon/2 FlexibleSUSY")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EFFHIGGSCOUPLINGS" << '\n'
            << FORMAT_RANK_THREE_TENSOR(25, 22, 22, (Abs(OBSERVABLES.eff_cp_higgs_photon_photon)), "Abs(effective H-Photon-Photon coupling)")
            << FORMAT_RANK_THREE_TENSOR(25, 21, 21, (Abs(OBSERVABLES.eff_cp_higgs_gluon_gluon)), "Abs(effective H-Gluon-Gluon coupling)")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices of
 * all given models in the SLHA object.
 *
 * @todo Use generic lambda instead of Set_spectrum in C++14
 *
 * @param models model classes
 */
template <class... Ts>
void cSMHdCKMRHNEFT_slha_io::set_spectrum(const std::tuple<Ts...>& models)
{
   Set_spectrum<cSMHdCKMRHNEFT_slha_io> ss(this);
   boost::fusion::for_each(models, ss);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void cSMHdCKMRHNEFT_slha_io::set_spectrum(const cSMHdCKMRHN<T>& model)
{
   set_spectrum(cSMHdCKMRHN_slha<cSMHdCKMRHN<T> >(model));
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::set_spectrum(const cSMHdCKMRHN_slha<Model>& model)
{
   const cSMHdCKMRHN_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void cSMHdCKMRHNEFT_slha_io::set_spectrum(const cSMHdCKM<T>& model)
{
   set_spectrum(cSMHdCKM_slha<cSMHdCKM<T> >(model));
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class Model>
void cSMHdCKMRHNEFT_slha_io::set_spectrum(const cSMHdCKM_slha<Model>& model)
{
   const cSMHdCKM_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);

   if (slha_io.get_modsel().quark_flavour_violated)
      set_ckm(model.get_ckm_matrix(), model.get_scale());

   if (slha_io.get_modsel().lepton_flavour_violated)
      set_pmns(model.get_pmns_matrix(), model.get_scale());
}

} // namespace flexiblesusy

#undef Pole
#undef PHYSICAL
#undef PHYSICAL_SLHA
#undef LOCALPHYSICAL
#undef MODEL
#undef MODELPARAMETER
#undef EXTRAPARAMETER
#undef OBSERVABLES
#undef LowEnergyConstant
#undef SCALES

#endif

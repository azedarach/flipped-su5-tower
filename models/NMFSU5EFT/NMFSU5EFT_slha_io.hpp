#ifndef NMFSU5EFT_SLHA_IO_H
#define NMFSU5EFT_SLHA_IO_H

#include "NMFSU5_mass_eigenstates.hpp"
#include "NMFSU5_model_slha.hpp"
#include "NMFSU5_info.hpp"
#include "NMFSU5_physical.hpp"

#include "cSMHdCKMRHNEFT_model_slha.hpp"

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

struct NMFSU5EFT_input_parameters;
class Spectrum_generator_settings;

template <class T>
class cSMHdCKM;

template <class T>
class cSMHdCKMRHN;

template <class T>
class NMFSU5;

struct NMFSU5EFT_scales {
   double HighScale{0.}, SUSYScale{0.}, LowScale{0.};
   double pole_mass_scale{0.};
};

class NMFSU5EFT_slha_io {
public:
   NMFSU5EFT_slha_io();

   void clear();

   void fill(softsusy::QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(NMFSU5EFT_input_parameters&) const;
   void fill(NMFSU5_mass_eigenstates&) const;
   template <class Model> void fill(NMFSU5_slha<Model>&) const;
   void fill(Physical_input&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void set_block(const std::string& str, SLHA_io::Position position = SLHA_io::back) { slha_io.set_block(str, position); }
   void set_blocks(const std::vector<std::string>& vec, SLHA_io::Position position = SLHA_io::back) { slha_io.set_blocks(vec, position); }
   template <class Model> void set_extra(const NMFSU5_slha<Model>&, const NMFSU5EFT_scales&, const NMFSU5_observables&);
   void set_input(const NMFSU5EFT_input_parameters&);
   void set_modsel(const SLHA_io::Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class... Ts> void set_spectrum(const std::tuple<Ts...>&);
   template <class Model> void set_spectrum(const cSMHdCKM_slha<Model>&);
   template <class T> void set_spectrum(const cSMHdCKM<T>&);
   template <class Model> void set_spectrum(const cSMHdCKMRHN_slha<Model>&);
   template <class T> void set_spectrum(const cSMHdCKMRHN<T>&);
   template <class Model> void set_spectrum(const NMFSU5_slha<Model>&);
   template <class T> void set_spectrum(const NMFSU5<T>&);
   void set_spectrum(const standard_model::Standard_model&);
   void set_spinfo(const Spectrum_generator_problems&);
   void set_spinfo(const Problems&);
   void set_spinfo(const std::vector<std::string>&, const std::vector<std::string>&);
   void set_print_imaginary_parts_of_majorana_mixings(bool);
   void write_to(const std::string&) const;
   void write_to_file(const std::string& file_name) const { slha_io.write_to_file(file_name); }
   void write_to_stream(std::ostream& ostr = std::cout) const { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(NMFSU5EFT_input_parameters&, int, double);
   static void fill_extpar_tuple(NMFSU5EFT_input_parameters&, int, double);
   static void fill_imminpar_tuple(NMFSU5EFT_input_parameters&, int, double);
   static void fill_imextpar_tuple(NMFSU5EFT_input_parameters&, int, double);

   template <class EFT, class IntermediateModel, class Model>
   static void fill_slhaea(SLHAea::Coll&, const cSMHdCKM_slha<EFT>&,
                           const cSMHdCKMRHN_slha<IntermediateModel>&,
                           const NMFSU5_slha<Model>&,
                           const softsusy::QedQcd&,
                           const NMFSU5EFT_input_parameters&,
                           const NMFSU5EFT_scales&, const NMFSU5_observables&);

   template <class EFT, class IntermediateModel, class Model>
   static SLHAea::Coll fill_slhaea(const cSMHdCKM_slha<EFT>&,
                                   const cSMHdCKMRHN_slha<IntermediateModel>&,
                                   const NMFSU5_slha<Model>&,
                                   const softsusy::QedQcd&,
                                   const NMFSU5EFT_input_parameters&,
                                   const NMFSU5EFT_scales&,
                                   const NMFSU5_observables&);

private:
   SLHA_io slha_io; ///< SLHA io class
   bool print_imaginary_parts_of_majorana_mixings;

   void set_extpar(const NMFSU5EFT_input_parameters&);
   void set_imminpar(const NMFSU5EFT_input_parameters&);
   void set_imextpar(const NMFSU5EFT_input_parameters&);
   void set_minpar(const NMFSU5EFT_input_parameters&);
   void set_mass(const cSMHdCKM_physical&, bool);
   void set_mass(const cSMHdCKMRHN_physical&, bool);
   void set_mass(const NMFSU5_physical&, bool);
   void set_mass(const standard_model::Standard_model_physical&);
   void set_mixing_matrices(const cSMHdCKM_physical&, bool);
   void set_mixing_matrices(const cSMHdCKMRHN_physical&, bool);
   void set_mixing_matrices(const NMFSU5_physical&, bool);
   void set_mixing_matrices(const standard_model::Standard_model_physical&);
   template <class Model> void set_model_parameters(const cSMHdCKM_slha<Model>&);
   template <class Model> void set_model_parameters(const cSMHdCKMRHN_slha<Model>&);
   template <class Model> void set_model_parameters(const NMFSU5_slha<Model>&);
   void set_model_parameters(const standard_model::Standard_model&);
   void set_ckm(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   void set_pmns(const Eigen::Matrix<std::complex<double>,3,3>&, double);
   double read_scale() const;
   void fill_drbar_parameters(NMFSU5_mass_eigenstates&) const;
   void fill_physical(NMFSU5_physical&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class Model>
void NMFSU5EFT_slha_io::fill(NMFSU5_slha<Model>& model) const
{
   fill(static_cast<NMFSU5_mass_eigenstates&>(model));
   fill_physical(model.get_physical_slha());
}

template <class EFT, class IntermediateModel, class Model>
void NMFSU5EFT_slha_io::fill_slhaea(
   SLHAea::Coll& slhaea, const cSMHdCKM_slha<EFT>& eft,
   const cSMHdCKMRHN_slha<IntermediateModel>& model,
   const NMFSU5_slha<Model>& high_scale_model,
   const softsusy::QedQcd& qedqcd, const NMFSU5EFT_input_parameters& input,
   const NMFSU5EFT_scales& scales, const NMFSU5_observables& observables)
{
   NMFSU5EFT_slha_io slha_io;
   const auto& problems = high_scale_model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_input(input);
   if (!error) {
      slha_io.set_spectrum(eft);
      slha_io.set_spectrum(model);
      slha_io.set_spectrum(high_scale_model);
      slha_io.set_extra(high_scale_model, scales, observables);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class EFT, class IntermediateModel, class Model>
SLHAea::Coll NMFSU5EFT_slha_io::fill_slhaea(
   const cSMHdCKM_slha<EFT>& eft,
   const cSMHdCKMRHN_slha<IntermediateModel>& model,
   const NMFSU5_slha<Model>& high_scale_model,
   const softsusy::QedQcd& qedqcd, const NMFSU5EFT_input_parameters& input,
   const NMFSU5EFT_scales& scales, const NMFSU5_observables& observables)
{
   SLHAea::Coll slhaea;
   NMFSU5EFT_slha_io::fill_slhaea(slhaea, eft, model, high_scale_model,
                                  qedqcd, input, scales, observables);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class Model>
void NMFSU5EFT_slha_io::set_model_parameters(const cSMHdCKM_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block EFT2gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("EFT2Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block_imag("IMEFT2Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("EFT2Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block_imag("IMEFT2Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("EFT2Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block_imag("IMEFT2Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block EFT2SM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Lambdax)), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("EFT2Kappa", ToMatrix(MODELPARAMETER(Kappa)), "Kappa", model.get_scale());
   slha_io.set_block_imag("IMEFT2Kappa", ToMatrix(MODELPARAMETER(Kappa)), "Kappa", model.get_scale());
   {
      std::ostringstream block;
      block << "Block EFT2HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
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
void NMFSU5EFT_slha_io::set_model_parameters(const cSMHdCKMRHN_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block EFT1gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "g1 * 0.7745966692414834")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("EFT1Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block_imag("IMEFT1Yu", ToMatrix(MODELPARAMETER(Yu_slha)), "Yu", model.get_scale());
   slha_io.set_block("EFT1Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block_imag("IMEFT1Yd", ToMatrix(MODELPARAMETER(Yd_slha)), "Yd", model.get_scale());
   slha_io.set_block("EFT1Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block_imag("IMEFT1Ye", ToMatrix(MODELPARAMETER(Ye_slha)), "Ye", model.get_scale());
   slha_io.set_block("EFT1Yv", MODELPARAMETER(Yv), "Yv", model.get_scale());
   slha_io.set_block_imag("IMEFT1Yv", MODELPARAMETER(Yv), "Yv", model.get_scale());
   slha_io.set_block("EFT1Mv", MODELPARAMETER(Mv), "Mv", model.get_scale());
   slha_io.set_block_imag("IMEFT1Mv", MODELPARAMETER(Mv), "Mv", model.get_scale());
   {
      std::ostringstream block;
      block << "Block EFT1SM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Lambdax)), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block EFT1HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
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
void NMFSU5EFT_slha_io::set_model_parameters(const NMFSU5_slha<Model>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(4, (MODELPARAMETER(gX)), "gX")
            << FORMAT_ELEMENT(5, (MODELPARAMETER(g5)), "g5")
         ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Y1", ToMatrix(MODELPARAMETER(Y1_slha)), "Y1", model.get_scale());
   slha_io.set_block_imag("ImY1", ToMatrix(MODELPARAMETER(Y1_slha)), "Y1", model.get_scale());
   slha_io.set_block("Y1Pr", ToMatrix(MODELPARAMETER(Y1Pr_slha)), "Y1Pr", model.get_scale());
   slha_io.set_block_imag("ImY1Pr", ToMatrix(MODELPARAMETER(Y1Pr_slha)), "Y1Pr", model.get_scale());
   slha_io.set_block("Y5b", ToMatrix(MODELPARAMETER(Y5b_slha)), "Y5b", model.get_scale());
   slha_io.set_block_imag("ImY5b", ToMatrix(MODELPARAMETER(Y5b_slha)), "Y5b", model.get_scale());
   slha_io.set_block("Y5bPr", ToMatrix(MODELPARAMETER(Y5bPr_slha)), "Y5bPr", model.get_scale());
   slha_io.set_block_imag("ImY5bPr", ToMatrix(MODELPARAMETER(Y5bPr_slha)), "Y5bPr", model.get_scale());
   slha_io.set_block("Y10", ToMatrix(MODELPARAMETER(Y10_slha)), "Y10", model.get_scale());
   slha_io.set_block_imag("ImY10", ToMatrix(MODELPARAMETER(Y10_slha)), "Y10", model.get_scale());
   slha_io.set_block("Y10Pr", ToMatrix(MODELPARAMETER(Y10Pr_slha)), "Y10Pr", model.get_scale());
   slha_io.set_block_imag("ImY10Pr", ToMatrix(MODELPARAMETER(Y10Pr_slha)), "Y10Pr", model.get_scale());

   {
      std::ostringstream block;
      block << "Block SU5LAMBDA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(lam1)), "lam1")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(lam2)), "lam2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(lam3)), "lam3")
            << FORMAT_ELEMENT(4, (MODELPARAMETER(lam3t)), "lam3t")
            << FORMAT_ELEMENT(5, (MODELPARAMETER(lam4)), "lam4")
            << FORMAT_ELEMENT(6, (MODELPARAMETER(lam4t)), "lam4t")
            << FORMAT_ELEMENT(7, (MODELPARAMETER(lam5)), "lam5")
            << FORMAT_ELEMENT(8, (MODELPARAMETER(lam5t)), "lam5t")
            << FORMAT_ELEMENT(9, (MODELPARAMETER(lam6)), "lam6")
            << FORMAT_ELEMENT(10, (MODELPARAMETER(lam6t)), "lam6t")
            << FORMAT_ELEMENT(11, Re(MODELPARAMETER(lam7)), "Re(lam7)")
            << FORMAT_ELEMENT(12, Re(MODELPARAMETER(lam8)), "Re(lam8)")
         ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block IMSU5LAMBDA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(11, Im(MODELPARAMETER(lam7)), "Im(lam7)")
            << FORMAT_ELEMENT(12, Im(MODELPARAMETER(lam8)), "Im(lam8)")
         ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block SU5ETA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, Re(MODELPARAMETER(eta1)), "Re(eta1)")
            << FORMAT_ELEMENT(2, Re(MODELPARAMETER(eta2)), "Re(eta2)")
            << FORMAT_ELEMENT(3, Re(MODELPARAMETER(eta3)), "Re(eta3)")
         ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block IMSU5ETA Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, Im(MODELPARAMETER(eta1)), "Im(eta1)")
            << FORMAT_ELEMENT(2, Im(MODELPARAMETER(eta2)), "Im(eta2)")
            << FORMAT_ELEMENT(3, Im(MODELPARAMETER(eta3)), "Im(eta3)")
         ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block SU5HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(m10sq)), "m10sq")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(m5sq)), "m5sq")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(m5Prsq)), "m5Prsq")
            << FORMAT_ELEMENT(4, Re(MODELPARAMETER(m12)), "Re(m12sq")
            << FORMAT_ELEMENT(5, Re(MODELPARAMETER(mu)), "Re(mu)")
            << FORMAT_ELEMENT(6, Re(MODELPARAMETER(muPr)), "Re(muPr)")
            << FORMAT_ELEMENT(7, (MODELPARAMETER(v)), "v")
            << FORMAT_ELEMENT(8, (MODELPARAMETER(vPr)), "vPr")
            << FORMAT_ELEMENT(9, (MODELPARAMETER(VG)), "VG")
         ;
      slha_io.set_block(block);
   }

   {
      std::ostringstream block;
      block << "Block IMSU5HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(4, Im(MODELPARAMETER(m12)), "Im(m12sq")
            << FORMAT_ELEMENT(5, Im(MODELPARAMETER(mu)), "Im(mu)")
            << FORMAT_ELEMENT(6, Im(MODELPARAMETER(muPr)), "Im(muPr)")
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
void NMFSU5EFT_slha_io::set_extra(
   const NMFSU5_slha<Model>& model, const NMFSU5EFT_scales& scales,
   const NMFSU5_observables& observables)
{

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
void NMFSU5EFT_slha_io::set_spectrum(const std::tuple<Ts...>& models)
{
   Set_spectrum<NMFSU5EFT_slha_io> ss(this);
   boost::fusion::for_each(models, ss);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void NMFSU5EFT_slha_io::set_spectrum(const NMFSU5<T>& model)
{
   set_spectrum(NMFSU5_slha<NMFSU5<T> >(model));
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class Model>
void NMFSU5EFT_slha_io::set_spectrum(const NMFSU5_slha<Model>& model)
{
   const NMFSU5_physical physical(model.get_physical_slha());
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
void NMFSU5EFT_slha_io::set_spectrum(const cSMHdCKMRHN<T>& model)
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
void NMFSU5EFT_slha_io::set_spectrum(const cSMHdCKMRHN_slha<Model>& model)
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
void NMFSU5EFT_slha_io::set_spectrum(const cSMHdCKM<T>& model)
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
void NMFSU5EFT_slha_io::set_spectrum(const cSMHdCKM_slha<Model>& model)
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

#include "NMFSU5EFT_slha_io.hpp"
#include "NMFSU5EFT_input_parameters.hpp"

#include "logger.hpp"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "config.h"

#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define PHYSICAL_SLHA(p) model.get_physical_slha().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) input.p
#define EXTRAPARAMETER(p) model.get_##p()
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;
#define LowEnergyConstant(p) Electroweak_constants::p

namespace flexiblesusy {

NMFSU5EFT_slha_io::NMFSU5EFT_slha_io()
   : slha_io()
   , print_imaginary_parts_of_majorana_mixings(false)
{
}

void NMFSU5EFT_slha_io::clear()
{
   slha_io.clear();
}

void NMFSU5EFT_slha_io::set_print_imaginary_parts_of_majorana_mixings(bool flag)
{
   print_imaginary_parts_of_majorana_mixings = flag;
}

/**
 * Stores the EXTPAR input parameters in the SLHA object.
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::set_extpar(const NMFSU5EFT_input_parameters& input)
{
   std::ostringstream extpar;

   extpar << "Block EXTPAR\n";
   extpar << FORMAT_ELEMENT(1, input.sign_delta_mAsq, "sign_delta_mAsq");
   extpar << FORMAT_ELEMENT(2, input.Lambda1IN, "Lambda1IN");
   extpar << FORMAT_ELEMENT(3, input.Lambda2IN, "Lambda2IN");
   extpar << FORMAT_ELEMENT(4, input.Lambda3IN, "Lambda3IN");
   extpar << FORMAT_ELEMENT(5, input.Lambda3tIN, "Lambda3tIN");
   extpar << FORMAT_ELEMENT(6, input.Lambda4IN, "Lambda4IN");
   extpar << FORMAT_ELEMENT(7, input.Lambda4tIN, "Lambda4tIN");
   extpar << FORMAT_ELEMENT(8, input.Lambda5IN, "Lambda5IN");
   extpar << FORMAT_ELEMENT(9, input.Lambda5tIN, "Lambda5tIN");
   extpar << FORMAT_ELEMENT(10, input.Lambda6IN, "Lambda6IN");
   extpar << FORMAT_ELEMENT(11, input.Lambda6tIN, "Lambda6tIN");
   extpar << FORMAT_ELEMENT(12, Re(input.Lambda7IN), "Re(Lambda7IN)");
   extpar << FORMAT_ELEMENT(13, Re(input.Lambda8IN), "Re(Lambda8IN)");
   extpar << FORMAT_ELEMENT(14, Re(input.Eta1IN), "Re(Eta1IN)");
   extpar << FORMAT_ELEMENT(15, Re(input.Eta2IN), "Re(Eta2IN)");
   extpar << FORMAT_ELEMENT(16, Re(input.Eta3IN), "Re(Eta3IN)");
   slha_io.set_block(extpar);
}

/**
 * Stores the IMMINPAR input parameters in the SLHA object.
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::set_imminpar(const NMFSU5EFT_input_parameters&)
{

}

/**
 * Stores the IMEXTPAR input parameters in the SLHA object.
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::set_imextpar(const NMFSU5EFT_input_parameters&)
{
   std::ostringstream extpar;

   extpar << "Block IMEXTPAR\n";
   extpar << FORMAT_ELEMENT(12, Im(input.Lambda7IN), "Im(Lambda7IN)");
   extpar << FORMAT_ELEMENT(13, Im(input.Lambda8IN), "Im(Lambda8IN)");
   extpar << FORMAT_ELEMENT(14, Im(input.Eta1IN), "Im(Eta1IN)");
   extpar << FORMAT_ELEMENT(15, Im(input.Eta2IN), "Im(Eta2IN)");
   extpar << FORMAT_ELEMENT(16, Im(input.Eta3IN), "Im(Eta3IN)");
   slha_io.set_block(extpar);
}

/**
 * Stores the MODSEL input parameters in the SLHA object.
 *
 * @param modsel struct of MODSEL parameters
 */
void NMFSU5EFT_slha_io::set_modsel(const SLHA_io::Modsel& modsel)
{
   slha_io.set_modsel(modsel);
}

/**
 * Stores the MINPAR input parameters in the SLHA object.
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::set_minpar(const NMFSU5EFT_input_parameters& input)
{
   std::ostringstream minpar;

   minpar << "Block MINPAR\n";

   slha_io.set_block(minpar);
}

/**
 * Stores all input parameters in the SLHA object.
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::set_input(const NMFSU5EFT_input_parameters& input)
{
   set_minpar(input);
   set_extpar(input);
   set_imminpar(input);
   set_imextpar(input);

   slha_io.set_block("UvIN", INPUTPARAMETER(UvInput), "UvInput");
   slha_io.set_block_imag("ImUvIN", INPUTPARAMETER(UvInput), "UvInput");
}

/**
 * Stores the additional physical input (FlexibleSUSYInput block) in
 * the SLHA object.
 *
 * @param input class of input
 */
void NMFSU5EFT_slha_io::set_physical_input(const Physical_input& input)
{
   slha_io.set_physical_input(input);
}

/**
 * Stores the settings (FlexibleSUSY block) in the SLHA object.
 *
 * @param settings class of settings
 */
void NMFSU5EFT_slha_io::set_settings(const Spectrum_generator_settings& settings)
{
   slha_io.set_settings(settings);
}

/**
 * Stores the SMINPUTS input parameters in the SLHA object.
 *
 * @param qedqcd class of Standard Model parameters
 */
void NMFSU5EFT_slha_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   slha_io.set_sminputs(qedqcd);
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void NMFSU5EFT_slha_io::set_spinfo(const Spectrum_generator_problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the spectrum generator information in the SPINFO block in
 * the SLHA object.
 *
 * @param problems struct with parameter point problems
 */
void NMFSU5EFT_slha_io::set_spinfo(const Problems& problems)
{
   set_spinfo(problems.get_problem_strings(), problems.get_warning_strings());
}

/**
 * Stores the given problems and warnings in the SPINFO block in the
 * SLHA object.
 *
 * @param problems vector of problem strings
 * @param warnings vector of warning strings
 */
void NMFSU5EFT_slha_io::set_spinfo(
   const std::vector<std::string>& problems,
   const std::vector<std::string>& warnings)
{
   std::ostringstream spinfo;
   spinfo << "Block SPINFO\n"
          << FORMAT_SPINFO(1, PKGNAME)
          << FORMAT_SPINFO(2, FLEXIBLESUSY_VERSION);

   for (const auto& s: warnings)
      spinfo << FORMAT_SPINFO(3, s);

   for (const auto& s: problems)
      spinfo << FORMAT_SPINFO(4, s);

   spinfo << FORMAT_SPINFO(5, NMFSU5_info::model_name)
          << FORMAT_SPINFO(9, SARAH_VERSION);

   slha_io.set_block(spinfo, SLHA_io::front);
}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void NMFSU5EFT_slha_io::set_mass(const cSMHdCKM_physical& physical,
                                 bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block EFT2MASS\n"
        << FORMAT_MASS(24, LOCALPHYSICAL(MVWp), "VWp")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(25, LOCALPHYSICAL(Mhh), "hh")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void NMFSU5EFT_slha_io::set_mass(const cSMHdCKMRHN_physical& physical,
                                 bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block EFT1MASS\n"
        << FORMAT_MASS(24, LOCALPHYSICAL(MVWp), "VWp")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(25, LOCALPHYSICAL(Mhh), "hh")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(8810012, LOCALPHYSICAL(MFv(3)), "Fv(4)")
         << FORMAT_MASS(8810014, LOCALPHYSICAL(MFv(4)), "Fv(5)")
         << FORMAT_MASS(8810016, LOCALPHYSICAL(MFv(5)), "Fv(6)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the particle masses in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_masses flag to indicate if Standard Model
 *    particle masses should be written as well
 */
void NMFSU5EFT_slha_io::set_mass(const NMFSU5_physical& physical,
                                 bool write_sm_masses)
{
   std::ostringstream mass;

   mass << "Block MASS\n"
        << FORMAT_MASS(42, LOCALPHYSICAL(MDelta(0)), "Delta(1)")
        << FORMAT_MASS(43, LOCALPHYSICAL(MDelta(1)), "Delta(2)")
        << FORMAT_MASS(44, LOCALPHYSICAL(MDelta(2)), "Delta(3)")
        << FORMAT_MASS(35, LOCALPHYSICAL(Mhh(1)), "hh(2)")
        << FORMAT_MASS(41, LOCALPHYSICAL(MVXY), "VXY")
        << FORMAT_MASS(32, LOCALPHYSICAL(MVZp), "VZp")
        << FORMAT_MASS(24, LOCALPHYSICAL(MVWp), "VWp")
   ;

   if (write_sm_masses) {
      mass
         << FORMAT_MASS(21, LOCALPHYSICAL(MVG), "VG")
         << FORMAT_MASS(25, LOCALPHYSICAL(Mhh(0)), "hh(1)")
         << FORMAT_MASS(1, LOCALPHYSICAL(MFd(0)), "Fd(1)")
         << FORMAT_MASS(3, LOCALPHYSICAL(MFd(1)), "Fd(2)")
         << FORMAT_MASS(5, LOCALPHYSICAL(MFd(2)), "Fd(3)")
         << FORMAT_MASS(2, LOCALPHYSICAL(MFu(0)), "Fu(1)")
         << FORMAT_MASS(4, LOCALPHYSICAL(MFu(1)), "Fu(2)")
         << FORMAT_MASS(6, LOCALPHYSICAL(MFu(2)), "Fu(3)")
         << FORMAT_MASS(11, LOCALPHYSICAL(MFe(0)), "Fe(1)")
         << FORMAT_MASS(13, LOCALPHYSICAL(MFe(1)), "Fe(2)")
         << FORMAT_MASS(15, LOCALPHYSICAL(MFe(2)), "Fe(3)")
         << FORMAT_MASS(12, LOCALPHYSICAL(MFv(0)), "Fv(1)")
         << FORMAT_MASS(14, LOCALPHYSICAL(MFv(1)), "Fv(2)")
         << FORMAT_MASS(16, LOCALPHYSICAL(MFv(2)), "Fv(3)")
         << FORMAT_MASS(8810012, LOCALPHYSICAL(MFv(3)), "Fv(4)")
         << FORMAT_MASS(8810014, LOCALPHYSICAL(MFv(4)), "Fv(5)")
         << FORMAT_MASS(8810016, LOCALPHYSICAL(MFv(5)), "Fv(6)")
         << FORMAT_MASS(22, LOCALPHYSICAL(MVP), "VP")
         << FORMAT_MASS(23, LOCALPHYSICAL(MVZ), "VZ")
      ;
   }

   slha_io.set_block(mass);

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void NMFSU5EFT_slha_io::set_mixing_matrices(const cSMHdCKM_physical& physical,
                                            bool write_sm_mixing_matrics)
{

   if (write_sm_mixing_matrics) {
      slha_io.set_block("EFT2UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("EFT2UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("EFT2UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("EFT2UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("EFT2UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("EFT2UERMIX", LOCALPHYSICAL(Ue), "Ue");
      slha_io.set_block("EFT2UVMIX", LOCALPHYSICAL(UV), "UV");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
      slha_io.set_block_imag("IMEFT2UVMIX", LOCALPHYSICAL(UV), "UV");
   }

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void NMFSU5EFT_slha_io::set_mixing_matrices(const cSMHdCKMRHN_physical& physical,
                                            bool write_sm_mixing_matrics)
{

   if (write_sm_mixing_matrics) {
      slha_io.set_block("EFT1UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("EFT1UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("EFT1UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("EFT1UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("EFT1UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("EFT1UERMIX", LOCALPHYSICAL(Ue), "Ue");
      slha_io.set_block("EFT1UVMIX", LOCALPHYSICAL(UV), "UV");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
      slha_io.set_block_imag("IMEFT1UVMIX", LOCALPHYSICAL(UV), "UV");
   }

}

/**
 * Stores the mixing matrices in the SLHA object.
 *
 * @param physical struct of physical parameters
 *
 * @param write_sm_mixing_matrics flag to indicate if Standard Model
 *    particle mixing matrices should be written as well
 */
void NMFSU5EFT_slha_io::set_mixing_matrices(const NMFSU5_physical& physical,
                                            bool write_sm_mixing_matrics)
{
   slha_io.set_block("UDELTAMIX", LOCALPHYSICAL(UDelta), "UDelta");
   slha_io.set_block("UH", LOCALPHYSICAL(UH), "UH");

   if (write_sm_mixing_matrics) {
      slha_io.set_block("EFT1UULMIX", LOCALPHYSICAL(Vu), "Vu");
      slha_io.set_block("EFT1UDLMIX", LOCALPHYSICAL(Vd), "Vd");
      slha_io.set_block("EFT1UURMIX", LOCALPHYSICAL(Uu), "Uu");
      slha_io.set_block("EFT1UDRMIX", LOCALPHYSICAL(Ud), "Ud");
      slha_io.set_block("EFT1UELMIX", LOCALPHYSICAL(Ve), "Ve");
      slha_io.set_block("EFT1UERMIX", LOCALPHYSICAL(Ue), "Ue");
      slha_io.set_block("EFT1UVMIX", LOCALPHYSICAL(UV), "UV");
   }

   if (print_imaginary_parts_of_majorana_mixings) {
      slha_io.set_block_imag("IMEFT1UVMIX", LOCALPHYSICAL(UV), "UV");
   }

}

void NMFSU5EFT_slha_io::set_ckm(
   const Eigen::Matrix<std::complex<double>,3,3>& ckm_matrix,
   double scale)
{
   slha_io.set_block("VCKM"  , ckm_matrix.real(), "Re(CKM)", scale);
   slha_io.set_block("IMVCKM", ckm_matrix.imag(), "Im(CKM)", scale);
}

void NMFSU5EFT_slha_io::set_pmns(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns_matrix,
   double scale)
{
   slha_io.set_block("VPMNS"  , pmns_matrix.real(), "Re(PMNS)", scale);
   slha_io.set_block("IMVPMNS", pmns_matrix.imag(), "Im(PMNS)", scale);
}

void NMFSU5EFT_slha_io::set_model_parameters(const standard_model::Standard_model& model)
{
   {
      std::ostringstream block;
      block << "Block SMGAUGE Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_g1() * standard_model_info::normalization_g1), "gY")
            << FORMAT_ELEMENT(2, (model.get_g2()), "g2")
            << FORMAT_ELEMENT(3, (model.get_g3()), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("SMYu", ToMatrix(model.get_Yu()), "Yu", model.get_scale());
   slha_io.set_block("SMYd", ToMatrix(model.get_Yd()), "Yd", model.get_scale());
   slha_io.set_block("SMYe", ToMatrix(model.get_Ye()), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SMSM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (model.get_mu2()), "mu2")
            << FORMAT_ELEMENT(2, (model.get_Lambdax()), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block SMHMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (model.get_v()), "v")
      ;
      slha_io.set_block(block);
   }
}

void NMFSU5EFT_slha_io::set_mass(const standard_model::Standard_model_physical& physical)
{
   std::ostringstream mass;

   mass << "Block SMMASS\n"
      << FORMAT_MASS(24, physical.MVWp, "VWp")
      << FORMAT_MASS(21, physical.MVG, "VG")
      << FORMAT_MASS(12, physical.MFv(0), "Fv(1)")
      << FORMAT_MASS(14, physical.MFv(1), "Fv(2)")
      << FORMAT_MASS(16, physical.MFv(2), "Fv(3)")
      << FORMAT_MASS(25, physical.Mhh, "hh")
      << FORMAT_MASS(1, physical.MFd(0), "Fd(1)")
      << FORMAT_MASS(3, physical.MFd(1), "Fd(2)")
      << FORMAT_MASS(5, physical.MFd(2), "Fd(3)")
      << FORMAT_MASS(2, physical.MFu(0), "Fu(1)")
      << FORMAT_MASS(4, physical.MFu(1), "Fu(2)")
      << FORMAT_MASS(6, physical.MFu(2), "Fu(3)")
      << FORMAT_MASS(11, physical.MFe(0), "Fe(1)")
      << FORMAT_MASS(13, physical.MFe(1), "Fe(2)")
      << FORMAT_MASS(15, physical.MFe(2), "Fe(3)")
      << FORMAT_MASS(22, physical.MVP, "VP")
      << FORMAT_MASS(23, physical.MVZ, "VZ")
      ;

   slha_io.set_block(mass);
}

void NMFSU5EFT_slha_io::set_mixing_matrices(const standard_model::Standard_model_physical& physical)
{
   slha_io.set_block("SMUULMIX", physical.Vu, "Vu");
   slha_io.set_block("SMUDLMIX", physical.Vd, "Vd");
   slha_io.set_block("SMUURMIX", physical.Uu, "Uu");
   slha_io.set_block("SMUDRMIX", physical.Ud, "Ud");
   slha_io.set_block("SMUELMIX", physical.Ve, "Ve");
   slha_io.set_block("SMUERMIX", physical.Ue, "Ue");
}

void NMFSU5EFT_slha_io::set_spectrum(const standard_model::Standard_model& model)
{
   const auto& physical = model.get_physical();

   set_model_parameters(model);
   set_mass(physical);
   set_mixing_matrices(physical);
}

/**
 * Write SLHA object to given output.  If output == "-", then the SLHA
 * object is written to std::cout.  Otherwise, output is interpreted
 * as a file name
 *
 * @param output "-" for cout, or file name
 */
void NMFSU5EFT_slha_io::write_to(const std::string& output) const
{
   if (output == "-")
      write_to_stream(std::cout);
   else
      write_to_file(output);
}

/**
 * Read (DR-bar) model parameter output scale from MODSEL entry 12
 */
double NMFSU5EFT_slha_io::get_parameter_output_scale() const
{
   return slha_io.get_modsel().parameter_output_scale;
}

/**
 * Read SLHA object from file
 *
 * @param file_name file name
 */
void NMFSU5EFT_slha_io::read_from_file(const std::string& file_name)
{
   slha_io.read_from_file(file_name);
}

/**
 * Read SLHA object from source
 *
 * calls SLHA_io::read_from_source()
 *
 * @param source source name
 */
void NMFSU5EFT_slha_io::read_from_source(const std::string& source)
{
   slha_io.read_from_source(source);
}

/**
 * Read SLHA object from stream
 *
 * @param istr stream name
 */
void NMFSU5EFT_slha_io::read_from_stream(std::istream& istr)
{
   slha_io.read_from_stream(istr);
}

/**
 * Fill struct of model input parameters from SLHA object (MINPAR,
 * EXTPAR and IMEXTPAR blocks)
 *
 * @param input struct of model input parameters
 */
void NMFSU5EFT_slha_io::fill(NMFSU5EFT_input_parameters& input) const
{
   SLHA_io::Tuple_processor minpar_processor = [&input, this] (int key, double value) {
      return fill_minpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor extpar_processor = [&input, this] (int key, double value) {
      return fill_extpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imminpar_processor = [&input, this] (int key, double value) {
      return fill_imminpar_tuple(input, key, value);
   };

   SLHA_io::Tuple_processor imextpar_processor = [&input, this] (int key, double value) {
      return fill_imextpar_tuple(input, key, value);
   };

   slha_io.read_block("MINPAR", minpar_processor);
   slha_io.read_block("EXTPAR", extpar_processor);
   slha_io.read_block("IMMINPAR", imminpar_processor);
   slha_io.read_block("IMEXTPAR", imextpar_processor);

   Eigen::Matrix<double, 3, 3> ReUvInput(Eigen::Matrix<double, 3, 3>::Zero());
   Eigen::Matrix<double, 3, 3> ImUvInput(Eigen::Matrix<double, 3, 3>::Zero());
   slha_io.read_block("UvIN", ReUvInput);
   slha_io.read_block("ImUvIN", ImUvInput);

   input.UvInput.real() = ReUvInput;
   input.UvInput.imag() = ImUvInput;
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 *
 * @param model model class to be filled
 */
void NMFSU5EFT_slha_io::fill_drbar_parameters(NMFSU5_mass_eigenstates& model) const
{
   model.set_gX(slha_io.read_entry("gauge", 4));
   model.set_g5(slha_io.read_entry("gauge", 5));

   {
      Eigen::Matrix<std::complex<double>,3,3> Y10;
      slha_io.read_block("Y10", Y10);
      model.set_Y10(Y10);
   }
   {
      Eigen::Matrix<std::complex<double>,3,3> Y10Pr;
      slha_io.read_block("Y10Pr", Y10Pr);
      model.set_Y10Pr(Y10Pr);
   }
   {
      Eigen::Matrix<std::complex<double>,3,3> Y5b;
      slha_io.read_block("Y5b", Y5b);
      model.set_Y5b(Y5b);
   }
   {
      Eigen::Matrix<std::complex<double>,3,3> Y5bPr;
      slha_io.read_block("Y5bPr", Y5bPr);
      model.set_Y5bPr(Y5bPr);
   }
   {
      Eigen::Matrix<std::complex<double>,3,3> Y1;
      slha_io.read_block("Y1", Y1);
      model.set_Y1(Y1);
   }
   {
      Eigen::Matrix<std::complex<double>,3,3> Y1Pr;
      slha_io.read_block("Y1Pr", Y1Pr);
      model.set_Y1Pr(Y1Pr);
   }

   model.set_lam1(slha_io.read_entry("SU5LAMBDA", 1));
   model.set_lam2(slha_io.read_entry("SU5LAMBDA", 2));
   model.set_lam3(slha_io.read_entry("SU5LAMBDA", 3));
   model.set_lam3t(slha_io.read_entry("SU5LAMBDA", 4));
   model.set_lam4(slha_io.read_entry("SU5LAMBDA", 5));
   model.set_lam4t(slha_io.read_entry("SU5LAMBDA", 6));
   model.set_lam5(slha_io.read_entry("SU5LAMBDA", 7));
   model.set_lam5t(slha_io.read_entry("SU5LAMBDA", 8));
   model.set_lam6(slha_io.read_entry("SU5LAMBDA", 9));
   model.set_lam6t(slha_io.read_entry("SU5LAMBDA", 10));
   model.set_lam7(std::complex<double>(
                     slha_io.read_entry("SU5LAMBDA", 11),
                     slha_io.read_entry("IMSU5LAMBDA", 11)));
   model.set_lam8(std::complex<double>(
                     slha_io.read_entry("SU5LAMBDA", 12),
                     slha_io.read_entry("IMSU5LAMBDA", 12)));
   model.set_eta1(std::complex<double>(
                     slha_io.read_entry("SU5ETA", 1),
                     slha_io.read_entry("IMSU5ETA", 1)));
   model.set_eta2(std::complex<double>(
                     slha_io.read_entry("SU5ETA", 2),
                     slha_io.read_entry("IMSU5ETA", 2)));
   model.set_eta3(std::complex<double>(
                     slha_io.read_entry("SU5ETA", 3),
                     slha_io.read_entry("IMSU5ETA", 3)));
   model.set_m10sq(slha_io.read_entry("SU5HMIX", 1));
   model.set_m5sq(slha_io.read_entry("SU5HMIX", 2));
   model.set_m5Prsq(slha_io.read_entry("SU5HMIX", 3));
   model.set_m12sq(std::complex<double>(
                      slha_io.read_entry("SU5HMIX", 4),
                      slha_io.read_entry("IMSU5HMIX", 4)));
   model.set_mu(std::complex<double>(
                   slha_io.read_entry("SU5HMIX", 5),
                   slha_io.read_entry("IMSU5HMIX", 5)));
   model.set_muPr(std::complex<double>(
                     slha_io.read_entry("SU5HMIX", 6),
                     slha_io.read_entry("IMSU5HMIX", 6)));
   model.set_v(slha_io.read_entry("SU5HMIX", 7));
   model.set_vPr(slha_io.read_entry("SU5HMIX", 8));
   model.set_VG(slha_io.read_entry("SU5HMIX", 9));


   model.set_scale(read_scale());
}

/**
 * Reads DR-bar parameters, pole masses and mixing matrices (in
 * Haber-Kane convention) from a SLHA output file.
 *
 * @param model model class to be filled
 */
void NMFSU5EFT_slha_io::fill(NMFSU5_mass_eigenstates& model) const
{
   fill_drbar_parameters(model);

   NMFSU5_physical physical_hk;
   fill_physical(physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

/**
 * Fill struct of extra physical input parameters from SLHA object
 * (FlexibleSUSYInput block)
 *
 * @param input struct of physical non-SLHA input parameters
 */
void NMFSU5EFT_slha_io::fill(Physical_input& input) const
{
   slha_io.fill(input);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings to be filled
 */
void NMFSU5EFT_slha_io::fill(Spectrum_generator_settings& settings) const
{
   slha_io.fill(settings);
}

void NMFSU5EFT_slha_io::fill_minpar_tuple(NMFSU5EFT_input_parameters& input,
                                          int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block MINPAR: " << key); break;
   }

}

void NMFSU5EFT_slha_io::fill_extpar_tuple(NMFSU5EFT_input_parameters& input,
                                          int key, double value)
{
   switch (key) {
   case 1: input.sign_delta_mAsq = Sign(value); break;
   case 2: input.Lambda1IN = value; break;
   case 3: input.Lambda2IN = value; break;
   case 4: inputl.Lambda3IN = value; break;
   case 5: input.Lambda3tIN = value; break;
   case 6: input.Lambda4IN = value; break;
   case 7: input.Lambda4tIN = value; break;
   case 8: input.Lambda5IN = value; break;
   case 9: input.Lambda5tIN = value; break;
   case 10: input.Lambda6IN = value; break;
   case 11: input.Lambda6tIN = value; break;
   case 12: input.Lambda7IN = std::complex<double>(value, Im(input.Lambda7IN)); break;
   case 13: input.Lambda8IN = std::complex<double>(value, Im(input.Lambda8IN)); break;
   case 14: input.Eta1IN = std::complex<double>(value, Im(input.Eta1IN)); break;
   case 15: input.Eta2IN = std::complex<double>(value, Im(input.Eta2IN)); break;
   case 16: input.Eta3IN = std::complex<double>(value, Im(input.Eta3IN)); break;
   default: WARNING("Unrecognized entry in block EXTPAR: " << key); break;
   }

}

void NMFSU5EFT_slha_io::fill_imminpar_tuple(NMFSU5EFT_input_parameters& /* input */,
                                            int key, double value)
{
   switch (key) {
   default: WARNING("Unrecognized entry in block IMMINPAR: " << key); break;
   }

}

void NMFSU5EFT_slha_io::fill_imextpar_tuple(NMFSU5EFT_input_parameters& input,
                                            int key, double value)
{
   switch (key) {
   case 12: input.Lambda7IN = std::complex<double>(Re(input.Lambda7IN), value); break;
   case 13: input.Lambda8IN = std::complex<double>(Re(input.Lambda8IN), value); break;
   case 14: input.Eta1IN = std::complex<double>(Re(input.Eta1IN), value); break;
   case 15: input.Eta2IN = std::complex<double>(Re(input.Eta2IN), value); break;
   case 16: input.Eta3IN = std::complex<double>(Re(input.Eta3IN), value); break;
   default: WARNING("Unrecognized entry in block IMEXTPAR: " << key); break;
   }

}

/**
 * Reads pole masses and mixing matrices from a SLHA output file to be filled.
 */
void NMFSU5EFT_slha_io::fill_physical(NMFSU5_physical& physical) const
{
   {
      DEFINE_PHYSICAL_PARAMETER(Vu);
      slha_io.read_block("UULMIX", Vu);
      LOCALPHYSICAL(Vu) = Vu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Vd);
      slha_io.read_block("UDLMIX", Vd);
      LOCALPHYSICAL(Vd) = Vd;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Uu);
      slha_io.read_block("UURMIX", Uu);
      LOCALPHYSICAL(Uu) = Uu;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ud);
      slha_io.read_block("UDRMIX", Ud);
      LOCALPHYSICAL(Ud) = Ud;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ve);
      slha_io.read_block("UELMIX", Ve);
      LOCALPHYSICAL(Ve) = Ve;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(Ue);
      slha_io.read_block("UERMIX", Ue);
      LOCALPHYSICAL(Ue) = Ue;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UV);
      slha_io.read_block("UVMIX", UV);
      LOCALPHYSICAL(UV) = UV;
   }

   LOCALPHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   LOCALPHYSICAL(Mhh) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   LOCALPHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   LOCALPHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   LOCALPHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   LOCALPHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   LOCALPHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   LOCALPHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   LOCALPHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   LOCALPHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   LOCALPHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   LOCALPHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   LOCALPHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   LOCALPHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   LOCALPHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   LOCALPHYSICAL(MFv)(3) = slha_io.read_entry("MASS", 8810012);
   LOCALPHYSICAL(MFv)(4) = slha_io.read_entry("MASS", 8810014);
   LOCALPHYSICAL(MFv)(5) = slha_io.read_entry("MASS", 8810016);
   LOCALPHYSICAL(MVWp) = slha_io.read_entry("MASS", 24);

}

/**
 * Reads the renormalization scales from all DR-bar parameter blocks.
 * If blocks with different scales are found the last scale is
 * returned and a warning is printed.
 *
 * @return common renormalization scale
 */
double NMFSU5EFT_slha_io::read_scale() const
{
   static const std::array<std::string, 8> drbar_blocks =
      { "gauge", "Y10", "Y10Pr", "Y5b", "Y5bPr", "Y1", "Y1Pr",
        "SU5LAMBDA", "SU5ETA", "SU5HMIX" };

   double scale = 0.;

   for (const auto& block: drbar_blocks) {
      const double block_scale = slha_io.read_scale(block);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

} // namespace flexiblesusy

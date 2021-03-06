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
 * @file cSMHdCKM_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

// File generated at Mon 5 Nov 2018 12:47:51

#ifndef cSMHdCKM_SLHA_H
#define cSMHdCKM_SLHA_H

#include "cSMHdCKM_input_parameters.hpp"
#include "cSMHdCKM_mass_eigenstates.hpp"
#include "cSMHdCKM_physical.hpp"

#include "ckm.hpp"
#include "linalg2.hpp"
#include "pmns.hpp"
#include "slha_io.hpp"
#include "wrappers.hpp"

#define LOCALPHYSICAL(p) physical.p
#define INPUTPARAMETER(p) this->input.p
#define MODELPARAMETER(p) this->p
#define PHYSICAL_SLHA(p) physical_slha.p
#define PHYSICAL_SLHA_REAL(p) Re(physical_slha.p)

namespace flexiblesusy {

/**
 * @class cSMHdCKM_slha<T>
 * @brief model class wrapper in SLHA convention
 *
 * @tparam Model model class to wrap
 */

template <class Model>
class cSMHdCKM_slha : public Model {
public:
   explicit cSMHdCKM_slha(const cSMHdCKM_input_parameters& input_ = cSMHdCKM_input_parameters());
   explicit cSMHdCKM_slha(const Model&, bool do_convert_masses_to_slha = true);
   cSMHdCKM_slha(const cSMHdCKM_slha&) = default;
   cSMHdCKM_slha(cSMHdCKM_slha&&) = default;
   virtual ~cSMHdCKM_slha() = default;
   cSMHdCKM_slha& operator=(const cSMHdCKM_slha&) = default;
   cSMHdCKM_slha& operator=(cSMHdCKM_slha&&) = default;

   virtual void clear() override;
   void convert_to_slha(); ///< converts pole masses and couplings to SLHA convention
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm_matrix() const { return ckm; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns_matrix() const { return pmns; }
   const cSMHdCKM_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   cSMHdCKM_physical& get_physical_slha(); ///< returns pole masses to SLHA convention
   void set_convert_masses_to_slha(bool); ///< allow/disallow for negative majoran fermion masses

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void print(std::ostream&) const override;

   double get_MVG_pole_slha() const { return PHYSICAL_SLHA(MVG); }
   double get_MHm_pole_slha() const { return PHYSICAL_SLHA(MHm); }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return PHYSICAL_SLHA(MFv); }
   double get_MFv_pole_slha(int i) const { return PHYSICAL_SLHA(MFv(i)); }
   double get_MAh_pole_slha() const { return PHYSICAL_SLHA(MAh); }
   double get_Mhh_pole_slha() const { return PHYSICAL_SLHA(Mhh); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return PHYSICAL_SLHA(MFd); }
   double get_MFd_pole_slha(int i) const { return PHYSICAL_SLHA(MFd(i)); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return PHYSICAL_SLHA(MFu); }
   double get_MFu_pole_slha(int i) const { return PHYSICAL_SLHA(MFu(i)); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return PHYSICAL_SLHA(MFe); }
   double get_MFe_pole_slha(int i) const { return PHYSICAL_SLHA(MFe(i)); }
   double get_MVWp_pole_slha() const { return PHYSICAL_SLHA(MVWp); }
   double get_MVP_pole_slha() const { return PHYSICAL_SLHA(MVP); }
   double get_MVZ_pole_slha() const { return PHYSICAL_SLHA(MVZ); }

   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd_pole_slha() const { return PHYSICAL_SLHA(Vd); }
   std::complex<double> get_Vd_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Vd(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud_pole_slha() const { return PHYSICAL_SLHA(Ud); }
   std::complex<double> get_Ud_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ud(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu_pole_slha() const { return PHYSICAL_SLHA(Vu); }
   std::complex<double> get_Vu_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Vu(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu_pole_slha() const { return PHYSICAL_SLHA(Uu); }
   std::complex<double> get_Uu_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Uu(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve_pole_slha() const { return PHYSICAL_SLHA(Ve); }
   std::complex<double> get_Ve_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ve(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue_pole_slha() const { return PHYSICAL_SLHA(Ue); }
   std::complex<double> get_Ue_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ue(i,k)); }
   const Eigen::Matrix<double,2,2>& get_ZZ_pole_slha() const { return PHYSICAL_SLHA(ZZ); }
   double get_ZZ_pole_slha(int i, int k) const { return PHYSICAL_SLHA(ZZ(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_UV_pole_slha() const { return PHYSICAL_SLHA(UV); }
   std::complex<double> get_UV_pole_slha(int i, int k) const { return PHYSICAL_SLHA(UV(i,k)); }

   const Eigen::Array<double,3,1>& get_Yu_slha() const { return Yu_slha; }
   double get_Yu_slha(int i) const { return Yu_slha(i); }
   const Eigen::Array<double,3,1>& get_Yd_slha() const { return Yd_slha; }
   double get_Yd_slha(int i) const { return Yd_slha(i); }
   const Eigen::Array<double,3,1>& get_Ye_slha() const { return Ye_slha; }
   double get_Ye_slha(int i) const { return Ye_slha(i); }
   const Eigen::Array<double,3,1>& get_Kappa_slha() const { return Kappa_slha; }
   double get_Kappa_slha(int i) const { return Kappa_slha(i); }



   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd_slha() const { return Vd_slha; }
   std::complex<double> get_Vd_slha(int i, int k) const { return Vd_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu_slha() const { return Vu_slha; }
   std::complex<double> get_Vu_slha(int i, int k) const { return Vu_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud_slha() const { return Ud_slha; }
   std::complex<double> get_Ud_slha(int i, int k) const { return Ud_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu_slha() const { return Uu_slha; }
   std::complex<double> get_Uu_slha(int i, int k) const { return Uu_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve_slha() const { return Ve_slha; }
   std::complex<double> get_Ve_slha(int i, int k) const { return Ve_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue_slha() const { return Ue_slha; }
   std::complex<double> get_Ue_slha(int i, int k) const { return Ue_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_UV_slha() const { return UV_slha; }
   std::complex<double> get_UV_slha(int i, int k) const { return UV_slha(i,k); }


private:
   cSMHdCKM_physical physical_slha{}; ///< contains the pole masses and mixings in slha convention
   Eigen::Matrix<std::complex<double>,3,3> ckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> pmns{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   bool convert_masses_to_slha{true};    ///< allow/disallow for negative majoran fermion masses
   Eigen::Array<double,3,1> Yu_slha{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Yd_slha{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Ye_slha{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Kappa_slha{Eigen::Array<double,3,1>::Zero()};

   Eigen::Matrix<std::complex<double>,3,3> Vd_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> UV_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};



   void calculate_ckm_matrix();
   void calculate_pmns_matrix();
   void convert_yukawa_couplings_to_slha();
   void convert_trilinear_couplings_to_slha();
   void convert_soft_squared_masses_to_slha();
};

template <class Model>
cSMHdCKM_slha<Model>::cSMHdCKM_slha(const cSMHdCKM_input_parameters& input_)
   : Model(input_)
{
}

/**
 * Copy constructor.  Copies from base class (model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 * @param do_convert_masses_to_slha whether to convert majorana
 *    fermion masses to SLHA convention (allow them to be negative)
 */
template <class Model>
cSMHdCKM_slha<Model>::cSMHdCKM_slha(const Model& model_,
                            bool do_convert_masses_to_slha)
   : Model(model_)
   , convert_masses_to_slha(do_convert_masses_to_slha)
{
   convert_to_slha();
}

template <class Model>
void cSMHdCKM_slha<Model>::clear()
{
   Model::clear();
   physical_slha.clear();
}

template <class Model>
void cSMHdCKM_slha<Model>::calculate_spectrum()
{
   Model::calculate_spectrum();
   convert_to_slha();
}

template <class Model>
void cSMHdCKM_slha<Model>::convert_to_slha()
{
   physical_slha = this->get_physical();

   if (convert_masses_to_slha)
      physical_slha.convert_to_slha();

   convert_yukawa_couplings_to_slha();
   calculate_ckm_matrix();
   calculate_pmns_matrix();
   convert_trilinear_couplings_to_slha();
   convert_soft_squared_masses_to_slha();
}

template <class Model>
void cSMHdCKM_slha<Model>::calculate_ckm_matrix()
{
   ckm = Vu_slha * Vd_slha.adjoint();
   CKM_parameters::to_pdg_convention(ckm, Vu_slha, Vd_slha, Uu_slha, Ud_slha);

}

template <class Model>
void cSMHdCKM_slha<Model>::calculate_pmns_matrix()
{
   const auto sign_delta_mAsq = INPUTPARAMETER(sign_delta_mAsq);

   Eigen::Matrix<double,3,3> neutrino_basis(Eigen::Matrix<double,3,3>::Zero());
   if (sign_delta_mAsq >= 0) {
      neutrino_basis = Eigen::Matrix<double,3,3>::Identity();
   } else {
      neutrino_basis(0,2) = 1.;
      neutrino_basis(1,0) = 1.;
      neutrino_basis(2,1) = 1.;
   }

   Eigen::Matrix<std::complex<double>,3,3> Uv(neutrino_basis.transpose() * UV_slha);
   pmns = Ve_slha * Uv.adjoint();
   PMNS_parameters::to_pdg_convention(pmns, Uv, Ve_slha, Ue_slha);
   UV_slha = neutrino_basis * Uv;
}

/**
 * Convert Yukawa couplings to SLHA convention
 */
template <class Model>
void cSMHdCKM_slha<Model>::convert_yukawa_couplings_to_slha()
{
   fs_svd(MODELPARAMETER(Yu), Yu_slha, Uu_slha, Vu_slha);
   fs_svd(MODELPARAMETER(Yd), Yd_slha, Ud_slha, Vd_slha);
   fs_svd(MODELPARAMETER(Ye), Ye_slha, Ue_slha, Ve_slha);
   fs_diagonalize_symmetric(this->get_Kappa(), Kappa_slha, UV_slha);

}

/**
 * Convert trilinear couplings to SLHA convention
 */
template <class Model>
void cSMHdCKM_slha<Model>::convert_trilinear_couplings_to_slha()
{

}

/**
 * Convert soft-breaking squared mass parameters to SLHA convention
 */
template <class Model>
void cSMHdCKM_slha<Model>::convert_soft_squared_masses_to_slha()
{

}

template <class Model>
const cSMHdCKM_physical& cSMHdCKM_slha<Model>::get_physical_slha() const
{
   return physical_slha;
}

template <class Model>
cSMHdCKM_physical& cSMHdCKM_slha<Model>::get_physical_slha()
{
   return physical_slha;
}

template <class Model>
void cSMHdCKM_slha<Model>::print(std::ostream& ostr) const
{
   Model::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
}

template <class Model>
void cSMHdCKM_slha<Model>::set_convert_masses_to_slha(bool flag)
{
   convert_masses_to_slha = flag;
}

template <class Model>
std::ostream& operator<<(std::ostream& ostr, const cSMHdCKM_slha<Model>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

#undef LOCALPHYSICAL
#undef INPUTPARAMETER
#undef MODELPARAMETER
#undef PHYSICAL_SLHA
#undef PHYSICAL_SLHA_REAL

#endif

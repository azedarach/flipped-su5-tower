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

#include "cSMHdCKMRHN_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

void cSMHdCKMRHN_physical::clear()
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

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void cSMHdCKMRHN_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MFv), LOCALPHYSICAL(UV));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void cSMHdCKMRHN_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MFv), LOCALPHYSICAL(UV));

}

Eigen::ArrayXd cSMHdCKMRHN_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

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

void cSMHdCKMRHN_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

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

Eigen::ArrayXd cSMHdCKMRHN_physical::get_masses() const
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

void cSMHdCKMRHN_physical::set_masses(const Eigen::ArrayXd& pars)
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

void cSMHdCKMRHN_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
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
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "UV = " << UV << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const cSMHdCKMRHN_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy

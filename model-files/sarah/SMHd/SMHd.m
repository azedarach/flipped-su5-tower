Off[General::spell];
Off[General::stop];

Model`Name = "SMHd";
Model`NameLaTeX = "SM with Y = -1/2 Higgs doublet";
Model`Authors = "D. Harries";
Model`Date = "2018-10-30";

Gauge[[1]] = {B,   U[1], hypercharge, g1, False};
Gauge[[2]] = {WB, SU[2],        left, g2, True};
Gauge[[3]] = {G,  SU[3],       color, g3, False};

FermionFields[[1]] = {q, 3, {uL, dL},  1/6, 2,  3};  
FermionFields[[2]] = {l, 3, {vL, eL}, -1/2, 2,  1};
FermionFields[[3]] = {d, 3, dLc,       1/3, 1, -3};
FermionFields[[4]] = {u, 3, uLc,      -2/3, 1, -3};
FermionFields[[5]] = {e, 3, eLc,         1, 1,  1};

ScalarFields[[1]] = {H,    1, {H0, Hm},  -1/2, 2, 1};

NameOfStates = {GaugeES, EWSB};

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC, {AddHC -> True}},
    {LagNoHC, {AddHC -> False}}
};

LagNoHC = -mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H
LagHC =  -(Yd H.d.q + Ye H.e.l + Yu u.q.conj[H] + WOp / 4 conj[H].l.conj[H].l);

DEFINITION[EWSB][GaugeSector] = {
    {{VB, VWB[3]}, {VP, VZ}, ZZ},
    {{VWB[1], VWB[2]}, {VWp, conj[VWp]}, ZW}
};

DEFINITION[EWSB][VEVs] = {
    {H0, {v, 1 / Sqrt[2]}, {Ah, \[ImaginaryI] / Sqrt[2]}, {hh, 1 / Sqrt[2]}}
};

DEFINITION[EWSB][MatterSector] = {
    {{{dL}, {dLc}}, {{DL, Vd}, {DR, Ud}}},
    {{{uL}, {uLc}}, {{UL, Vu}, {UR, Uu}}},
    {{{eL}, {eLc}}, {{EL, Ve}, {ER, Ue}}},
    {{vL}, {VL, UV}}
};

DEFINITION[GaugeES][DiracSpinors] = {
    Fd1 -> {dL, 0},
    Fd2 -> {0, dLc},
    Fu1 -> {uL, 0},
    Fu2 -> {0, uLc},
    Fe1 -> {eL, 0},
    Fe2 -> {0, eLc},
    Fv1 -> {vL, 0}
};

DEFINITION[EWSB][DiracSpinors] = {
    Fd -> {DL, conj[DR]},
    Fe -> {EL, conj[ER]},
    Fu -> {UL, conj[UR]},
    Fv -> {VL, 0}
};

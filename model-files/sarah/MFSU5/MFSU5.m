Off[General::spell];
Off[General::stop];

Model`Name = "MFSU5";
Model`NameLateX = "Minimal flipped $SU(5)$";
Model`Authors = "D. Harries";
Model`Date = "2018-10-09";

Gauge[[1]] = {BX, U[1], u1x, gX, False};
Gauge[[2]] = {A5, SU[5], five, g5, False};

ScalarFields[[1]] = {Phi10, 1, phi10, QXPhi10, 10};
ScalarFields[[2]] = {Phi5, 1, phi5, QXPhi5, 5};

NameOfStates = { GaugeES };

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC,   { AddHC -> True }},
    {LagNoHC, { AddHC -> False }}
};

LagHC = -(mu5 Phi10.Phi10.Phi5 / 8);

LagNoHC = -(m10sq conj[Phi10].Phi10 / 2
            + m5sq conj[Phi5].Phi5
            + lam1 conj[Phi10].Phi10.conj[Phi10].Phi10 / 4
            + lam2 conj[Phi10].Phi10.conj[Phi10].Phi10 / 4
            + lam3 conj[Phi5].Phi5.conj[Phi5].Phi5
            + lam4 conj[Phi5].Phi5.conj[Phi10].Phi10 / 2
            + lam5 conj[Phi5].Phi10.conj[Phi10].Phi5);

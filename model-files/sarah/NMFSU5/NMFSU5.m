Off[General::spell];
Off[General::stop];

Model`Name = "NMFSU5";
Model`NameLaTeX = "Next-to-minimal flipped $SU(5)$";
Model`Authors = "D. Harries";
Model`Date = "2018-10-02";

Gauge[[1]] = {BX, U[1], u1x, gX, False};
Gauge[[2]] = {A5, SU[5], five, g5, False};

ScalarFields[[1]] = {Phi10, 1, phi10, QXPhi10, 10};
ScalarFields[[2]] = {Phi5, 1, phi5, QXPhi5, 5};
ScalarFields[[3]] = {Phi5p, 1, phi5p, QXPhi5p, 5};

NameOfStates = {GaugeES};

DEFINITION[GaugeES][LagrangianInput] = {
(*    {LagHC,   {AddHC -> True}}, *)
    {LagNoHC, {AddHC -> False}}
};
(*
LagHC = -(m12sq conj[Phi5].Phi5p
          + mu5 Phi10.Phi10.Phi5 / 8 + mu5p Phi10.Phi10.Phi5p / 8
          + eta1 conj[Phi5].Phi5.conj[Phi5].Phi5p
          + eta2 conj[Phi5].Phi5p.conj[Phi5].Phi5p
          + eta3 conj[Phi5].Phi5p.conj[Phi5p].Phi5p
          + lam7 conj[Phi5].Phi5p.conj[Phi10].Phi10 / 2
          + lam8 conj[Phi5].Phi10.conj[Phi10].Phi5p);
*)

LagNoHC = -(m10sq conj[Phi10].Phi10 / 2
            + m5sq conj[Phi5].Phi5 + m5psq conj[Phi5p].Phi5p
            + lam1 conj[Phi10].Phi10.conj[Phi10].Phi10 / 4
            + lam2 conj[Phi10].Phi10.conj[Phi10].Phi10 / 4
            + lam3 conj[Phi5].Phi5.conj[Phi5].Phi5
            + lam3t conj[Phi5p].Phi5p.conj[Phi5p].Phi5p
            + lam6 conj[Phi5].Phi5p.conj[Phi5p].Phi5
            + lam6t conj[Phi5].Phi5.conj[Phi5p].Phi5p
            + lam4 conj[Phi5].Phi5.conj[Phi10].Phi10 / 2
            + lam4t conj[Phi5p].Phi5p.conj[Phi10].Phi10 / 2
            + lam5 conj[Phi5].Phi10.conj[Phi10].Phi5
            + lam5t conj[Phi5p].Phi10.conj[Phi10].Phi5p);

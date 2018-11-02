ParameterDefinitions = {
    {g1,        { Description -> "Hypercharge-Coupling"}},
    {g2,        { Description -> "Left-Coupling"}},
    {g3,        { Description -> "Strong-Coupling"}},
    {AlphaS,    {Description -> "Alpha Strong"}},
    {e,         { Description -> "electric charge"}},

    {Gf,        { Description -> "Fermi's constant"}},
    {aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

    {Yu,        { Description -> "Up-Yukawa-Coupling",
                  DependenceNum -> None}},
    {Yd,        { Description -> "Down-Yukawa-Coupling",
                  DependenceNum -> None}},
    {Ye,        { Description -> "Lepton-Yukawa-Coupling",
                  DependenceNum -> None}},
    {Yv,        { Description -> "Neutrino-Yukawa-Coupling",
                  Real -> False,
                  LaTeX -> "Y_{\\nu}",
                  LesHouches -> Yv,
                  OutputName -> Yv,
                  DependenceNum -> None}},

    {Mv,       { Description -> "Right-handed neutrino Majorana mass",
                 Real -> False,
                 LaTeX -> "M_{\\nu}",
                 LesHouches -> Mv,
                 OutputName -> Mv,
                 DependenceNum -> None}},

    {mu2,        { Description -> "SM Mu Parameter",
                   OutputName->m2SM}},
    {\[Lambda],  { Description -> "SM Higgs Selfcouplings",
                   DependenceNum -> Mass[hh]^2/(v^2)}},
    {v,          { Description -> "EW-VEV",
                   DependenceNum -> None,
                   DependenceSPheno -> None,
                   OutputName -> vvSM}},
    {mH2,        { Description -> "SM Higgs Mass Parameter"}},

    {ThetaW,    { Description -> "Weinberg-Angle",
                  DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]  }},

    {ZZ, {Description -> "Photon-Z Mixing Matrix"}},
    {ZW, {Description -> "W Mixing Matrix",
           Dependence ->   1/Sqrt[2] {{1, 1},
                      {\[ImaginaryI],-\[ImaginaryI]}} }},


    {Vu,        {Description -> "Left-Up-Mixing-Matrix"}},
    {Vd,        {Description -> "Left-Down-Mixing-Matrix"}},
    {Uu,        {Description -> "Right-Up-Mixing-Matrix"}},
    {Ud,        {Description -> "Right-Down-Mixing-Matrix"}},
    {Ve,        {Description -> "Left-Lepton-Mixing-Matrix"}},
    {Ue,        {Description -> "Right-Lepton-Mixing-Matrix"}},
    {UV,        {Description -> "Neutrino-Mixing-Matrix"}}
};

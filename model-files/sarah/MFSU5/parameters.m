ParameterDefinitions = {
    {g5, { Description -> "SU(5) gauge coupling",
           Real -> True,
           LesHouches -> {GAUGE, 4},
           LaTeX -> "g_5",
           OutputName -> g5
         } },
    {gX,   { Description -> "U(1)_X gauge coupling",
             Real -> True,
             LesHouches -> {GAUGE, 6},
             LaTeX -> "g_X",
             OutputName -> gX
           } },
    {lam1, { Description -> "(10-10)-(10-10) quartic coupling",
             Real -> True,
             LesHouches -> {SU5MIX, 1},
             LaTeX -> "\\lambda_1",
             OutputName -> lam1
           } },
    {lam2, { Description -> "(10-10-10-10) quartic coupling",
             Real -> True,
             LesHouches -> {SU5MIX, 2},
             LaTeX -> "\\lambda_2",
             OutputName -> lam2
           } },
    {lam3, { Description -> "(5-5-5-5) quartic coupling",
             Real -> True,
             LesHouches -> {SU5MIX, 3},
             LaTeX -> "\\lambda_3",
             OutputName -> lam3
           } },
    {lam4, { Description -> "(10-10)-(5-5) quartic coupling",
             Real -> True,
             LesHouches -> {SU5MIX, 5},
             LaTeX -> "\\lambda_4",
             OutputName -> lam4
           } },
    {lam5, { Description -> "(5-10-10-5) quartic coupling",
             Real -> True,
             LesHouches -> {SU5MIX, 7},
             LaTeX -> "\\lambda_5",
             OutputName -> lam5
           } },
    {m10sq,  { Description -> "Scalar 10 squared mass",
             Real -> True,
             LesHouches -> {SU5MIX, 16},
             LaTeX -> "m_{10}^2",
             OutputName -> m10sq
           } },
    {m5sq,  { Description -> "Scalar 5 squared mass",
             Real -> True,
             LesHouches -> {SU5MIX, 17},
             LaTeX -> "m_{5}^2",
             OutputName -> m5sq
           } },
    {mu5,  { Description -> "Scalar 10-10-5 mixing",
             Real -> False,
             LesHouches -> {SU5MIX, 20},
             LaTeX -> "\\mu",
             OutputName -> mu5
           } },
    {QXPhi10, { Description -> "10-plet scalar U(1)_X charge",
                Real -> True,
                LesHouches -> {U1XCHARGE, 1},
                LaTeX -> "Q^X_{10_S}",
                OutputName -> QXPhi10
              } },
    {QXPhi5, { Description -> "5-plet scalar U(1)_X charge",
                Real -> True,
                LesHouches -> {U1XCHARGE, 2},
                LaTeX -> "Q^X_{5_S}",
                OutputName -> QXPhi5
              } }
};

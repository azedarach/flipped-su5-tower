set title "cSMHdCKMRHN renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='cSMHdCKMRHN_rgflow.dat'

plot for [i=2:96+1] filename using 1:(column(i)) title columnhead(i)

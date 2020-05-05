set term pdfcairo
set outp "Polymers_Gyration.pdf"

set title "Chain Gyrationradius"

set xlabel "Stiffness"
set ylabel "Gyrationradius"
plot "Average-chain.txt" using 1 : 2 w lp title "rg", "Average-chain.txt" using 1 : 3 w lp title "rg_{pos}"
plot "Average-chain.txt" using 1 : ($2 / $3)  w lp title "rg / rg_{pos}"
plot "Average-chain.txt" using 1 : 3  w lp title "rg_{pos}"


set title "Chain Asphericity"
set xlabel "Stiffness"
set ylabel "Asphericity"
plot "Average-chain.txt" using 1 : 4 w lp title "b", "Average-chain.txt" using 1 : (($4) / ($2)**2) w lp title "b/Rg^2", "Average-chain.txt" using 1 : (($4) / ($3)**2) w lp title "b/Rg_{pos}^2"

set title "Chain Prolateness"
set xlabel "Stiffness"
set ylabel "Prolateness"
plot "Average-chain.txt" using 1 : 5 w lp title "S*"




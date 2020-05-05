set term pdfcairo
set outp "Polymers_Gyration.pdf"


#Chain
set title "Chain Gyrationradius"

set xlabel "Stiffness"
set ylabel "Gyrationradius"
plot "Average-chain.txt" using 1 : 2 w lp title "rg", "Average-chain.txt" using 1 : 3 w lp title "rg_{pos}"


set title "Chain Asphericity"
set xlabel "Stiffness"
set ylabel "Asphericity"
plot "Average-chain.txt" using 1 : 4 w lp title "b"
plot "Average-chain.txt" using 1 : (($4) / ($2)**2) w lp title "b/Rg^2", "Average-chain.txt" using 1 : (($4) / ($3)**2) w lp title "b/Rg_{pos}^2"

set title "Chain Prolateness"
set xlabel "Stiffness"
set ylabel "Prolateness"
plot "Average-chain.txt" using 1 : 5 w lp title "S*"



#Ring
set title "Ring Gyrationradius"

set xlabel "Stiffness"
set ylabel "Gyrationradius"
plot "Average-ring.txt" using 1 : 2 w lp title "rg", "Average-ring.txt" using 1 : 3 w lp title "rg_{pos}"


set title "Ring Asphericity"
set xlabel "Stiffness"
set ylabel "Asphericity"
plot "Average-ring.txt" using 1 : 4 w lp title "b"
plot "Average-ring.txt" using 1 : (($4) / ($2)**2) w lp title "b/Rg^2", "Average-ring.txt" using 1 : (($4) / ($3)**2) w lp title "b/Rg_{pos}^2"

set title "Ring Prolateness"
set xlabel "Stiffness"
set ylabel "Prolateness"
plot "Average-ring.txt" using 1 : 5 w lp title "S*"


#Crossovers


set title "Gyrationradii"
set xlabel "Stiffness"
set ylabel "Gyrationradius"
plot "Average-chain.txt" using 1 : 2 w lp title "chain", "Average-ring.txt" using 1 : 2 w lp title "ring","Average-star.txt" using 1 : 2 w lp title "star"


set title "Asphericities/Rg^2"
set xlabel "Stiffness"
set ylabel "Asphericity"
plot "Average-chain.txt" using 1 : (($4) / ($2)**2) w lp title "chain", "Average-ring.txt" using 1 : (($4) / ($2)**2) w lp title "ring","Average-star.txt" using 1 : (($4) / ($2)**2) w lp title "star"


set title "Prolatenesses"
set xlabel "Stiffness"
set ylabel "Prolateness"
plot "Average-chain.txt" using 1 : 5 w lp title "chain", "Average-ring.txt" using 1 : 5 w lp title "ring","Average-star.txt" using 1 : 5 w lp title "star"


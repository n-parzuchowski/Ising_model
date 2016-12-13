#!/usr/bin/env python

from matplotlib.pylab import *

fx = open( 'several_trials.dat' ,'r' )

temp  = []
trials =[[],[],[],[],[],[]] 
for line in fx:
  broken = line.strip().split()
  temp.append(float(broken[0]))

  for indx in range(1,7):
    trials[indx-1].append(float(broken[indx]))

for indx in range(1,7):
  plot(temp,trials[indx-1])

xlabel("Temperature")
ylabel("Magnetization")

text(0.05,0.8,"ISING MODEL FERROMAGNETISM",size=10)
text(0.05,0.7,"$H = -J \sum_{<ij>} \sigma_i \sigma_j$",size=10)
text(0.05,0.6,"$J = 0.15$, initial configuration: Ferromagnet",size=10)
fx.close()
savefig("Curie_Temp.pdf")


show()

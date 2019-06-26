# CardioFAN
DOI: 10.5281/zenodo.1807129
https://doi.org/10.5281/zenodo.1807129

Accompannied article to cite:
Seyed Vahedein, Y., Liberson, A.S., 2019. CardioFAN: open source platform for noninvasive assessment of pulse transit time and pulsatile flow in hyperelastic vascular networks. Biomech Model Mechanobiol. https://doi.org/10.1007/s10237-019-01163-z

***Finite volume monotone code for blood flow hemodynamics simulation and "Pulse Transit Time" (PTT) calculation!

## Please Read the following instructions:
-In "main.m", choose between TVD LAX_Wendroff and regular LAX_Wendroff algorithms by **KTVD** parameter to 1 or 0, respectively. 
 RUN_NET_LW is a LAX_Wendroff based code
 RUN_NET_TVD is a TVD LAX_Wendroff based code

-Viscohyperelastic version of the code use the designated folder and execute RUN_NET_VISCOHYPER.

-The code is ready to run for 55 vessels (parameters were set using data from Sherwin and Alastruey et al. 2012), 37 vessels (Matthys et 
 al 2007) and 26 vessels cases (Alastruey et al. 2016). For **arbitrary geometry** and properties please follow the instructions below.

-User can switch between these by setting **NVESSEL** to any of these options (55,37,26)! for arbitrary case the following functions 
 should be altered.

-**INPUT.m** is used to define all the areterial model properties, such as, number of vessels, number of mesh elements, and the vessel 
 mechanical properties. All the instructions are written as comments in front of each parameter assigned inside INPUT function.

-**LaNrCrC.m** is where you define all the properties related to the network, such as, NODE-CONNECT MATRIX (connecting vessels by nodes), 
 Vessel Lengths, Vessel Radii, Resistance, Compliances and etc.

-**FLOWINPUT.m** is where you can define your own flow input or change the currently available inputs (Alastruey and Figuera et al.)

-**UFUN.m** is to define velocity profile assignment method. Ex. use new velocity profile or assign the fitted values given in FLOWINPUT (options are ready to use)

-**PTTcalc.m** calculates PTT starting from the desired cardiac cycle. It will calculate PTT considering if the PTT calculation would 
 start at any point during the cycle and gives a plot of all PTTs as a function of starting point and a plot of flow waveform at the 
 desired cycle.  
  If user encounters errors running the code, it is due to small total TIME at INPUT and can be fixed by one of the following methods:
  1) Adjust the number of time-steps in line 45 of **PTTcalc.m** function. Ex. change it from 5000 to smaller value! 
  2) Increase the total calculation TIME at **INPUT.m**
  3) Turn off PTTcalculation by putting it to 0 at **INPUT.m**

<G.U.I WILL BE ADDED SOON>

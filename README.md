# CardioFAN
Finite volume-based monotone code for blood flow hemodynamics simulations!

Choose between TVD LAX_Wendroff and regular LAX_Wendroff codes by selecting KTVD=1 or 0 in "main.m".
For Viscohyperelastic go to the designated folder and execute RUN_NET_VISCOHYPER.

The code is ready to run for 55 vessels (Sherwina and Alastruey et al. 2012), 37 vessels (Matthys et al 2007) and 26 vessels cases (Alastruey et al. 2016)

User can check these by setting NVESSEL to any of these options (55,37,26)! for arbitrary case the following functions should be altered.

RUN_NET_LW is a LAX_Wendroff based code
RUN_NET_TVD is a TVD LAX_Wendroff based code

INPUT is the file that the user needs to define all the properties, such as, number of vessels, number of mesh elements, and the vessel properties. All the instructions are written as comments in front of the parameters assigned in INPUT function.

LaNrCrC is where you define all the properties related to the network, such as, NODE-CONNECT MATRIX (connecting vessels by nodes), Vessel Lengths, Vessel Radii, Resistance, Compliances and etc.

FLOWINPUT is where you can define your own flow input or change the currently available inputs picked from Alastruey and Figuera et al.

UFUN is the place to define velocity profile or assign the fitted values given in FLOWINPUT (options are ready to use)

PTTcalc will calculate PTT starting from the desired cardiac cycle. It will calculate PTT considering if the PTT calculation would start at any point during the cycle and gives a plot of all PTTs as a function of starting point and a plot of flow waveform at the desired cycle.  
If user encounters errors running the code, it is due to small total TIME at INPUT and can be fixed by one of the following methods:
1) Adjust the number of time-steps in line 45 of PTTcalc function. Ex. change it from 5000 to smaller value! 
2) Increase the total calculation TIME at INPUT
3) Turn off PTTcalculation by putting it to 0 at INPUT


<G.U.I WILL BE ADDED SOON>

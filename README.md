# CardioFAN
Finite volume-based monotone code for blood flow hemodynamics simulations!

Choose between TVD LAX_Wendroff and regular LAX_Wendroff codes by selecting KTVD=1 or 0 in "main.m".

The code is ready to run for 55 vessels (Sherwina and Alastruey et al. 2012), 37 vessels (Matthys et al 2007) and 26 vessels cases (Alastruey et al. 2016)

User can check these by setting NVESSEL to any of these options (55,37,26)! for arbitrary case the following functions should be altered.

RUN_NET_LW is a LAX_Wendroff based code
RUN_NET_TVD is a TVD LAX_Wendroff based code

INPUT is the file that the user needs to define all the properties, such as, number of vessels, number of mesh elements, and the vessel properties.

LaNrCrC is where you define all the properties related to the network, such as, NODE-CONNECT MATRIX (connecting vessels by nodes), Vessel Lengths, Vessel Radii, Resistance, Compliances and etc.

FLOWINPUT is where you can define your own flow input or change the currently available inputs picked from Alastruey and Figuera et al.

PTTcalc will calculate PTT in from a starting time step (which varries from the initial starttime (defined at INPUT) plus 5000 time steps). 
If user encounters errors running the code, it is due to small total TIME at INPUT and can be fixed by one of the following methods:
1) Decrease the number of time-steps from 5000 in PTT 
2) Increase the total calculation TIME at INPUT
3) Turn off PTTcalculation by putting it to 0 at INPUT

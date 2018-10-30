function  MAIN
%% put KTVD=0 if you would like to run LAX_WENDROFF code for cardiovascular simulation
%% put KTVD=1 if you would like to run TVD LAX_WENDROFF code for cardiovascular simulation

%% select TVD (1) or LW (0)
%======================================
KTVD=0;                  %default is 0


%% Auto Starting the LW or TVD code!
'Dont forget to edit INPUT function for arbitrary network or for changing the number of vessels'

if KTVD==1
  RUN_NET_TVD
else 
  RUN_NET_LW
end

end
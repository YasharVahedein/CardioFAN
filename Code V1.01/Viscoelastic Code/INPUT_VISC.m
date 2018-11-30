function [SCALAR,ARRAY1D,PTTfunc,AP,AU,AREA,AREAZ,CMK, XX,NODE_CONNECT,Pfit]=INPUT_VISC(KTVD)
%**************
% MIT License
% 
% Copyright (c) 2018 <Yashar Seyed Vahedein, Alexander Liberson>
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% The user is recommended to reference the first released publication based on this code:
% 
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%**************

%================================================================
% PARAMETER DESCRIPTION
%================================================================
%SCALAR - array of all scalar inputs
%ARRAY1D - compound array of all 1D input arrays
%AP,AU - 2D array, AP(1:NVESSEL, 1:NCELLMAX) - cell centered
%NODE_CONNECT - node connectivity matrix,
%NODE_CONNECT=[ALLNODES(1:NNODE),ivup,ivdn1,ivdn2], ivup, ivdn - global numbers
%for upstream and doownstrean vessel segment, ALLNODES - Node number
% IF USING PTT CALCULATOR, MAKE SURE TO DEFINE THE START AND FINISH TIME OF CALCULATION
% OTHERWISE YOU WILL RECEIVE AN ERROR

%% DEFINE and/or Check all the Input Properties
%=========================================================
%*Model properties (SHOULD BE DEFINED BY USER)
%* calculation time + COURANT NUMBER + MESHING + Numerical Scheme
TIME=2;                %Total Time of Calculation %s
T=0.25;                % Necessary to change for arbitrary inlet sin function
dtau=1;                %time-step size initialize
CFL=0.2;              %CFL number, for TVD use lower CFL (use 0.4)
NcellT=450;            %NUMBER OF ALL CELLS (for TVD use more mesh elements)
Nmin=6;                %Min # of cells for any vessel %Nmin=8 for TVD
NVESSEL=37;            % ***** NUMBER OF VESSELS ****** %

%* For initialization(DEFINED BASED ON EACH CASE)
INITopt=0;             % initialize with Linear Acoustic Model
if KTVD==1
    INITopt=1;         % If using TVD automatically initializes 
end
Pzero=0;               % Extramural pressure
Pd=0;                  % Diastlic pressure added to initial pressure option 2

%* PULSE TRANSIT TIME (PTT) calculation on arbitrary path
PTTcalculation=0;      % 1 for turning on PTT calculation, 0 turn off
if PTTcalculation==1   % if 1 then define the properties of arbitrary path and start time along the waveform
    STARTtime=10000;   % Global start time for PTT calculation (caution on large values)
    FINISHtime=17000;  % Global finish time (caution on large values)
    pttSTAvessel1=1;   % Starting vessel number
    pttFINvessel2=3;  % Final vessel number
    LLLmax=0.285;      %for calculated value based on vessel number use 0
end

%* Visual settings: Geometry Schematics Line Thickness Thresholds + Plot Acoustic Solution Results
threshold=0.007;       % m
threshold1=0.005;      % m
acousticplot=0;        % Acoustic model PLOTs Yes=1, No=0

%%* If not using 26,55 or 37 vessels cases, uncomment & define lines 65-79
%     a=0;                     %*****FUNG parameter*****% a=0:means regular Hook & Laplace Law
%     Period=0.96;             % Duration/Period of heart beat ic output  %s   -  ALASTUEY 2016
%     Rho=1060;                % 26 Vessel density of the flow %kg/m^3
%     mju=3.5e-3;              % Dynamic Viscosity   %Pa*s
%     r0=(12.4e-3);            % Root Radius %m -  ALASTUEY 2016
%     Pout=9200;               % Pressure fitted    %Pa
%     Rcorrect=0.6;            %outlet resistance correction for 26 vessels case Alastruey-2016
%     nonconstant=0;           % 1 for nonconstant and 0 for constant CMK
%   %* Terminal boundary model parameters
%     RCR=1;                   %1 for RCR and 0 for R model
%     Rcorrect=1.00001;        %55 and 37 vessel cases -outlet resistance correction for 26 vessels case Alastruey-2016
%   %*Tune the following exit R and C multipliers if needed with method given in Xiao 2014 - Int.J.Numer.Meth.Biomed.Engng
%     Rmult=1;                 %resistance multiplier , set 1 for no multiplier
%     Cmult=1;                 %compliance multiplier , set 1 for no multiplier
%     correctcompliance=1;     %correcting compliance using total resistance
% %   correctcompliance=0;       %turn off correcting compliance effect

%* 26, 55 and 37 Case specific properties
Rcorrect=1.00001;        %55 and 37 vessel cases -outlet resistance correction
if NVESSEL==26
    BestP=0;           % For best PRESSURE approx put to 1 - for best AREA approx put to zero
    a=8;               % %*****FUNG parameter*****%
    Pout=9200;         % Pressure fitted (Best Pressure Alastruey 2016)  
    Pzero=10000;
    TPeriod=0.98;       % Duration/Period of heart beat ic output  %s   -  ALASTUEY 2016
    Rho=1060;          % 26 Vessel density of the flow %kg/m^3
    mju=3.5e-3;        % Dynamic Viscosity   %Pa*s
    r0=(12.4e-3);      % Root Radius %m -  ALASTUEY 2016
    Rcorrect=0.6;      %outlet resistance correction for 26 vessels case Alastruey-2016
    Rmult=1;                 %resistance multiplier , set 1 for no multiplier
    Cmult=1.5;               %compliance multiplier , set 1 for no multiplier
    nonconstant=0;     % 1 for nonconstant and 0 for constant CMK
    correctcompliance=1;     %correcting compliance using total resistance
    RCR=1;
    if BestP==0
        Pout=4400;    % cappilary pressure (Best Area Alastruey 2016)
        a=2;          % Fung parameter 2
        nonconstant=0;
        Rmult=3.3;  % FOR cappilary - resistance multiplier
        Cmult=1.6;  % FOR cappilary - compliance multiplier
    end
    
elseif NVESSEL==55
    a=1.3;
    NcellT=400;
    TPeriod=0.95;         % 55 vessel
    Rho=1050;             % 55 vessel density
    mju=4e-3;             % Dynamic Viscosity   %Pa*s
    r0=(15.4e-3);         % Root Radius %m
    Pout=0;               % Cappilary pressure, or 1330
    Pzero=10000;
    Pd=1330;              % ambiant P0 pressure 
    nonconstant=1;
    RCR=1;
    Rmult=1;              %resistance multiplier , set 1 for no multiplier
    Cmult=1;              %compliance multiplier , set 1 for no multiplier
    correctcompliance=0;  %Best case 0
    
elseif NVESSEL==37
    a=0;              %%*****FUNG parameter*****%
    NcellT=300;
    TPeriod=0.82;     % 37 vessel
    Rho=1050;         % 37 vessel density
    mju=2.5e-3;       % Dynamic Viscosity   %Pa*s
    r0=(14.3e-3);     % Root Radius %m
    Pout=426.56;      % 0
    Pzero=8800;
    nonconstant=1;
    RCR=1;            %1 for RCR and 0 for R model
    Rmult=1;                 %resistance multiplier , set 1 for no multiplier
    Cmult=1;                 %compliance multiplier , set 1 for no multiplier
    correctcompliance=1;
end

%*Calculating Viscous Flow Parameters
nju=mju/Rho;      %Kinematic Viscosity %Pa*s*m^3/Kg
KR=-22*pi*nju;    %Coefficient for the viscous source term
% KR=0;           % Uncomment for Inviscid Flow




%% DEFINE Function related to Flow(t)
%=========================================================
A0=pi*r0^2;       %Root Reference Area   %m^2
[A,AU1,Pfit]=FLOWINPUT(TIME,TPeriod,A0,NVESSEL);


%% DEFINE Lengths, Angles, NODE_CONNECTivity matrix
%% DEFINE Terminal Resistance and Compliance Arrays
%% DEFINE Radii, Area, CMK, Thickness
%=========================================================
[LL,ALF,NODE_CONNECT,ARZ1,ARZ4,AREAZ1,AREAZ4,CMK1,CMK4,h,OMEGA,RESIST,COMPLIANCE]=LaNrCrC_VISC(NVESSEL,Pout,nonconstant);
%LL  - Vessel lengths
%ALF - Vessel angles of orientations with respect to coordinate system
%NODE_CONNECT(inode=1:NNODE,ivup,ivdn1,ivdn2)
'size(NODE_CONNECT) from INPUT','size(NODE_CONNECT)'

% maximum length for PTT
if PTTcalculation==1 && LLLmax==0 
    LLLmax=sum(pttSTAvessel1:pttSTAvessel2);
end

%% Automatically finds inlet and exit vessels
%=========================================================
ALLNODES=NODE_CONNECT(:,1);          %all nodes set
NNODE=length(ALLNODES);              %total number of all nodes
j=0;                                 %value set to zero
for i=1:NNODE                        %total nodes numbers
    ii=ALLNODES(i);                  %physical node number
    if(NODE_CONNECT(ii,2)==0)        %no upstream vessel
        IVinlet=NODE_CONNECT(ii,3);  %inlet vessel number found !
    end
    if(NODE_CONNECT(ii,3)==0 && NODE_CONNECT(ii,4)==0)
        j=j+1;
        vesselnum=NODE_CONNECT(ii,2);
        IVexit(j)=NODE_CONNECT(ii,2);%Exit vessel number found
    end
end
IVinlet, IVexit                      %Print Exit and Inlet Vessel Numbers

% Check of Exit Vessels
%========================
NEXIT=length(IVexit);    %number of exit vessels
length(RESIST)          %number of Resistances defined
length(IVexit)          %number of exit BCs


%% Automatically plots network schematics
%================================================================
% PLOT Arterial GEOMETRY
GEOMETRY(LL,ALF,NVESSEL,NODE_CONNECT,ARZ1,ARZ4,NNODE,threshold,threshold1)


%% Automatically MESHING, Time-Step Size and Vessel Geometry Interpolation
%========================================================================
%Distribute cells across each vessel
Ltotal=sum(LL);   %total length of the network

for iv=1:NVESSEL  %uniform distribution of cells across each vessel
    NCELL(iv)=ceil(NcellT/Ltotal*LL(iv));
    if(NCELL(iv)<Nmin)
        NCELL(iv)=Nmin;
    end
    HMESH(iv)=LL(iv)/NCELL(iv);
    Yold=[AREAZ1(iv)  AREAZ4(iv)];
    XXold=[0 LL(iv)];
    XXnew=linspace(0,LL(iv),NCELL(iv));
    a0=interp1(XXold,Yold,XXnew,'linear');
    AREAZ(iv, 1:NCELL(iv))=a0;
    Yold=[CMK1(iv)  CMK4(iv)];
    CMK(iv,1:NCELL(iv))=interp1(XXold,Yold,XXnew,'linear');
    %=========================================================
    % PLOT VARIABLE Cross Section
    %=========================================================
    %     figure
    %     plot(XXnew, AREAZ(iv,1:NCELL(iv)));
    %     hold on
    %=========================================================
    % In case of discontineous properties
    %=========================================================
    %     CMK(iv,1:NCELL(iv)/2)=5.21; Think about it
    %     CMK(iv,NCELL(iv)/2:NCELL(iv))=5;
    %=========================================================
    % Calculate Young Modulus for viscoelastic case
    %=========================================================
%         E(iv,1:NCELL(iv))=(3*Rho*sqrt(AREAZ(iv,1:NCELL(iv))/pi).*CMK(iv,1:NCELL(iv)).^2)/(4*h(iv));
%         OMEGA(iv,1:NCELL(iv))=(phi(iv)/(Rho^2)).*(E(iv,1:NCELL(iv)))./((CMK(iv,1:NCELL(iv)).^2).*sqrt(AREAZ(iv,1:NCELL(iv))));

    %=========================================================
    % In Case of Adding Bending Moments Problem
    %=========================================================
    %           %% Bending Moments Calcs
    %           EB(iv,1)=(2*E(iv,1)-E(iv,2));
    %           D(iv,1)=(EB(iv,1)*h(iv)^3)/(12*(1-0.5^2));
    %           for I=1:NCELL(iv)-1
    %              EB(iv,I+1)=(E(iv,I)+E(iv,I+1))/2;
    %              D(iv,I+1)=(EB(iv,I+1)*h(iv)^3)/(12*(1-0.5^2));
    %           end
    %           EB(iv,NCELL(iv)+1)=(2*E(iv,NCELL(iv))-E(iv,NCELL(iv)-1));
    %           D(iv,NCELL(iv)+1)=(EB(iv,NCELL(iv)+1)*h(iv)^3)/(12*(1-0.5^2));
    
    %=========================================================
    % Time-step size recalculation
    %=========================================================
    dtau=min(dtau,CFL*min(HMESH(iv)./CMK(iv,:)));
end

%=========================================================
% Number of Time-Steps
%=========================================================
NTAU= TIME/dtau;



%% Initialization
%========================================================================
% Step 1
%========================================================================
NCELLmax=max(NCELL);  %need to specify the length of 2D array
XX(1:NVESSEL, 1:NCELLmax)=0;    %centered coordinates
%GRID=========================================
for iv=1:NVESSEL
    XX1=linspace(HMESH(iv)/2,LL(iv)-HMESH(iv)/2,NCELL(iv));
    XX(iv,1:NCELL(iv))=XX1;
end



%========================================================================
% option 1 - Zero U, P and AREA=Pre-stress AREA
%========================================================================
AU=zeros(NVESSEL,NCELLmax);
AP=zeros(NVESSEL,NCELLmax)+Pzero;
AREA(1:NVESSEL,1:NCELLmax)=AREAZ(1:NVESSEL,1:NCELLmax);

%========================================================================
% COMPLIANCE CALCULATION - SPECIFIC FOR ALASTRUEY 2016
%========================================================================
[COMPLIANCEFINAL,RC,RP,RT]=RCRparameters(NVESSEL,LL,CMK,AREAZ,XX,NCELL,IVexit,NEXIT,RESIST,COMPLIANCE,Rho,Rmult,Cmult,correctcompliance,Rcorrect);


if INITopt==1
    %======================================================================
    % option 2 - Initialization with accoustics solution
    %======================================================================
    %     close all;
    % Prepare to send for acoustic solution
    SCALAR0=[Rho, KR, NVESSEL, NNODE,KTVD,CFL,dtau,NTAU/(TIME/TPeriod)*4,IVinlet,NEXIT,T,TPeriod,TIME,RCR,Pzero];
    ARRAY1D0=[NCELL, HMESH, LL, ALF, IVexit, RT,RC,RP,COMPLIANCEFINAL];
    % Calculate acoustic solution
    [p(1:NVESSEL, 1:NCELLmax),u(1:NVESSEL, 1:NCELLmax)]=ACOU_INIT_VISC(SCALAR0,ARRAY1D0,AP,AU,AREAZ,CMK,XX,NODE_CONNECT,Pfit,Rmult,Cmult,Pout,acousticplot); %%%%%%IMPORTANT
    AU=u(1:NVESSEL, 1:NCELLmax);
    AP=p(1:NVESSEL, 1:NCELLmax)+Pd;
    for iv=1:NVESSEL
        %         AREA(iv, 1:NCELLmax)=AREAZ(iv,1:NCELLmax)*(AP(iv,1:NCELLmax)/(2*Rho*CMK(iv,1:NCELLmax).^2)+1)^2;
        AREA(iv, 1:NCELLmax)=AREAZ(iv,1:NCELLmax)*(((AP(iv,1:NCELLmax)-Pzero)/(2*Rho*CMK(iv,1:NCELLmax).^2))*(1-a*((AP(iv,1:NCELLmax)-Pzero)/(2*Rho*(CMK(iv,1:NCELLmax).^2)))^2)+1).^2;
    end
end

%% output
%======================================================================
% INPUT DATA FOR THE NONLINEAR CODE
%======================================================================
'SCALAR from INPUT'
SCALAR=[Rho,KR, NVESSEL, NNODE,KTVD,CFL,dtau,NTAU,IVinlet,NEXIT,T,TPeriod,TIME,a,Pout,RCR,Pzero];
'ARRAY1D  from INPUT'
ARRAY1D=[NCELL, HMESH, LL, ALF, OMEGA, IVexit, RT,RC,RP,COMPLIANCEFINAL];
if PTTcalculation==1
    PTTfunc=[PTTcalculation,STARTtime,FINISHtime,pttSTAvessel1,pttFINvessel2,LLLmax];
else
    PTTfunc=[PTTcalculation];
end
end






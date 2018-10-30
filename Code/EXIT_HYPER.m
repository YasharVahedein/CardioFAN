function [AUBN,AUBNp, APBN,APBNp, AREABN,AREABNp, FBN, HBN, SBN]=EXIT_HYPER(APN,APBNp, AUN,AUBNp,A0N,AREAN,AREABNp, rho, dtau,t1, l, cmk, KR,RT,RC,RP,CT,a,Pout,RCR) %add AU1, dtau, l, AREA1,AUNp, AREANp,APNp
%**************
% Copyright (c) Rochester Institute of Technology (RIT) 2018 <Yashar Seyed Vahedein, Alexander Liberson>
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
%RC,RP,R(or RESIST) - Conduit, Peripheral and Total Resistances of arteries
%CC,CP,C(or COMPLIANCE) - Conduit, Peripheral and Total Compliance of arteries
%CT total corrected compliance
%AREABN, AREABNp - Area at the exit edge in current and previous time-step
%APBN, APBNp - Pressure at the exit edge in current and previous time-step
%AUBN, AUBNP - Velocity at the exit edge in current and previous time-step
%YN Forward moving linear characteristic wave
%alpha, beta and gamma - substitutions for terms in RCR formula
%FBN,SBN - Flux and Source terms at exit boundary
%HBN - Jacobian at the exit boundary
if RCR==0
    
    %% Matthys - 37 vessels nonlinear
    %         %% LINEAR - WORKS BEST
    RES=(RT-rho*cmk/(A0N))/(RT+rho*cmk/(A0N));
    %         RES=0;                                                   %non-reflective
    YN=AUN+APN/rho/cmk;
    APBN=(1+RES)*YN/2*rho*cmk;
    AUBN=(1-RES)*YN/2;
    AREABN=A0N*((APBN/(2*rho*(cmk^2))*(1-a*(APBN/(2*rho*(cmk^2)))^2))+1)^2;
    
elseif RCR==1
    %% If Resistances not available
    %    phi=2^(-1/3);                                  %IF R  not available
    %    lambda=1/(2^phi^2);                            %IF R  not available
    %    R=RESIST; %(lambda/((2*phi^4)-lambda))*       %IF R  not available
    %    C0=4.12*10^(-12);                              %IF CT not available
    %    Ct=((2*lambda*phi^3)/(1-(2*lambda*phi^3)))*C0; %IF CT not available
    %    C=COMPLIANCE; % use next from Matthys 1e-11; %
    %    L0= (RESIST^2)*C0;                             %IF RCLR method used
    %    L=(lambda/((2*phi^2)-lambda))*L0;              %IF RCLR method used
    
    % Calculations
    alpha=RC/dtau+(RC/RP+1)/(CT);
    beta=1/(dtau)+1/(RP*CT);
    gamma=AUBNp*AREABNp*RC/dtau-APBNp/(dtau)-Pout/(RP*CT);
    YN=AUN+APN/rho/cmk;
    AUBN=(rho*cmk*beta*YN+gamma)/(alpha*AREABNp+beta*rho*cmk);
    APBN=(YN-AUBN)*rho*cmk;
    AREABN=A0N*((APBN/(2*rho*(cmk^2))*(1-a*(APBN/(2*rho*(cmk^2)))^2))+1)^2; %AREABN=A0N*(APBN/(2*rho*cmk^2)+1)^2;
    AUBNp=AUBN;
    AREABNp=AREABN;
    APBNp=APBN;
end

%% FLUX and Jacobians
FBN=[AUBN*AREABN; AUBN^2/2+APBN/rho];
HBN=[AUBN AREABN; (cmk^2)/sqrt(A0N*AREABN)*exp(a*(sqrt(AREABN/A0N)-1)^2)*(1+2*a*(sqrt(AREABN/A0N)-1)^2) AUBN];
SBN=[0; KR*AUBN/AREABN];
end


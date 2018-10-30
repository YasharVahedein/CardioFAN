function [ AUB1, APB1, AREAB1, FB1, HB1, SB1] = INLET_HYPER(AU1, AREA1,A01, rho,dtau, t,TT, cmk,KR,Pfit,a)
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
%AREABN, AREABNp - Area at the exit edge in current and previous time-step
%APBN, APBNp - Pressure at the exit edge in current and previous time-step
%AUB1, APB1, AREAB1 - Velocity,pressure and area at the edge of the inlet vessel
%alpha - calculated for hyperelastic boundary wave propagation
%FB1,SB1 - Flux and Source terms at the edge inlet boundary
%HB1 - Jacobian at the exit boundary

    AUB1=UFUN_NONLIN(Pfit,t,TT,AREA1);  
    alpha=sqrt(AREA1*(cmk^2/sqrt(A01*AREA1))*exp(a*(sqrt(AREA1/A01)-1)^2)*(1+2*a*(sqrt(AREA1/A01)-1)^2));  %%define a
    AREAB1=(AUB1-AU1)/(alpha/AREA1)+AREA1;
    APB1=2*rho*cmk^2*(sqrt(AREAB1/A01)-1)*exp(a*(sqrt(AREAB1/A01)-1)^2);  
    FB1=[AUB1*AREAB1; AUB1^2/2+APB1/rho];
    HB1=[AUB1 AREAB1; (cmk^2)/sqrt(A01*AREAB1)*exp(a*(sqrt(AREAB1/A01)-1)^2)*(1+2*a*(sqrt(AREAB1/A01)-1)^2) AUB1];
    SB1=[0; KR*AUB1/AREAB1];
end

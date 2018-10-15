function [APB1N,AUB1N,AREAB1N,APB21,AUB21,AREAB21,APB31,AUB31,AREAB31,FB1N,FB21,FB31,HB1N,HB21,HB31,SB1N,SB21,SB31]=...
    BIFURC_HYPER(P1N,U1N,P21,U21,P31,U31,roc1,roc2,roc3,A1,A01N,A2,A021,A3,A031,cmk1N,cmk21,cmk31,kr1,kr2,kr3,rho,a)
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
%Calculates Riemann variables at bifurcations
%P1N,U1N - the last cell variables in the vessel 1
%P21,U21 - the first cell variables at the vessel 2
%P31,U31 - the first cell variables at the vessel 3
%PP1,UU1 - Riemann variables at the exit of the vessel 1
%PP2,UU2 - Riemann variables at the inlet to the vessel 2
%PP3,UU3 - Riemann variables at the inlet to the vessel 3
%if ibifurk=3 - fork; ibifurk=2 - area change only

I1=U1N+P1N/roc1;
 I2=U21-P21/roc2;
  I3=U31-P31/roc3;

PP=(I1*A1-I2*A2-I3*A3)/(A1/roc1+A2/roc2+A3/roc3);
APB1N=PP; APB21=PP; APB31=PP;

% AREAB1N=A01N*(APB1N/(2*rho*cmk1N^2)+1)^2; 
AREAB1N=A01N*((APB1N/(2*rho*(cmk1N^2))*(1-a*(APB1N/(2*rho*(cmk1N^2)))^2))+1)^2;
%  AREAB21=A021*(APB21/(2*rho*cmk21^2)+1)^2; 
 AREAB21=A021*((APB21/(2*rho*(cmk21^2))*(1-a*(APB21/(2*rho*(cmk21^2)))^2))+1)^2;
%   AREAB31=A031*(APB31/(2*rho*cmk31^2)+1)^2; 
  AREAB31=A031*((APB31/(2*rho*(cmk31^2))*(1-a*(APB31/(2*rho*(cmk31^2)))^2))+1)^2;

AUB1N=I1-APB1N/roc1;
 AUB21=I2+APB21/roc2;
  AUB31=I3+APB31/roc3;

FB1N=[AUB1N*AREAB1N; AUB1N^2/2+APB1N/rho];
HB1N=[AUB1N AREAB1N; (cmk1N^2)/sqrt(A01N*AREAB1N)*exp(a*(sqrt(AREAB1N/A01N)-1)^2)*(1+2*a*(sqrt(AREAB1N/A01N)-1)^2) AUB1N];
SB1N=[0; kr1*AUB1N/AREAB1N];
    FB21=[AUB21*AREAB21; AUB21^2/2+APB21/rho];
    HB21=[AUB21 AREAB21; (cmk21^2)/sqrt(A021*AREAB21)*exp(a*(sqrt(AREAB21/A021)-1)^2)*(1+2*a*(sqrt(AREAB21/A021)-1)^2) AUB21];
    SB21=[0; kr2*AUB21/AREAB21];
        FB31=[AUB31*AREAB31; AUB31^2/2+APB31/rho];
        HB31=[AUB31 AREAB31; (cmk31^2)/sqrt(A031*AREAB31)*exp(a*(sqrt(AREAB31/A031)-1)^2)*(1+2*a*(sqrt(AREAB31/A031)-1)^2) AUB31];
        SB31=[0; kr3*AUB31/AREAB31];
end


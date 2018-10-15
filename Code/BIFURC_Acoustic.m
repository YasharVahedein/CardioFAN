function [PP1,UU1,PP2,UU2,PP3,UU3]=...
    BIFURC_Acoustic(P1N,U1N,P21,U21,P31,U31,roc1,roc2,roc3,A1,A2,A3)
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
PP1=PP; PP2=PP; PP3=PP;
UU1=I1-PP1/roc1;
UU2=I2+PP2/roc2;
UU3=I3+PP3/roc3;
end


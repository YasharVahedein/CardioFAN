function [ PP, UU ] = INLET_Acoustic(P1, U1, A1, rhoc, t ,T,TT,Pfit)
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

%PRESSURE SPECIFIED AT THE INLET
%===========================================
% %PP,UU - Riemann variables at the inlet
% %P1,U1 - cell centered variables at the first cell
% %t - current node time +half a time step
% PP=PFUN(t);
% UU=U1+(PP-P1)/rhoc;

%INLET IS FIXED
%==================
% PP=P1;
% UU=U1;

%VELOCITY SPECIFIED AT THE INLET
 
UU=UFUN_NONLIN(Pfit,t,TT,A1);
PP=(UU-U1)*rhoc+P1;

% 'General: from Inlet: t,Tperiod,UU,U1,rho,cmk,P1'
% t,TT,UU,U1,rhoc,P1
return

end

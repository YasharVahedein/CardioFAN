function U=UFUN_NONLIN(Pfit,t,TT,A1)
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
%=================================
%Inlet velocity 
%t  - current time
%T  - length of active Cardiac Output
%TT - period of Cardiac output
%U  - Velocity at the boundary

if(t>TT), t=rem(t,TT); end

% rem(x,y) is x - n.*y where n = fix(x./y) if y ~= 0. 
% rem(23,5)=3

%% option 1
% % R0=(15.4e-3);        %m
% % A0=pi*R0^2;       %m^2
% QMAX=500; %mL/s; Flow rate
% QMAX=QMAX*0.001;  % L/s
% QMAX=QMAX*0.001;  % m^3/s
% 
% UMAX=QMAX/A1;

% tau=0.35;  arc=0.3; %0.25; 
% if t<arc
%     U=UMAX*sin(pi*t/arc)^1.3;
% elseif t<tau
%     U=UMAX/12*(-sin(6.8*pi*(t-arc)/tau)^(3/4));
% else
%     U=0;
% end

%% option2 - Matthys variable CS
%     U=polyval(Pfit,t).*(0.000001/(A1)); %7.4506e-04, 5.8965e-04

%% option3 - Find velocity from Fitted Flow(t)
    U=Pfit(t).*(0.000001/(A1)); %Q changed to U

end

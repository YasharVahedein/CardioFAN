function [AFMh]=TVD_Nonlinear(AU,AREA,A0,ncell,rho,dtau,dx,cmk,R,Lambda, KTVD) %AREAB1,AUB1,AREABN,AUBN,
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
%N number of cells in the current vessel
%Q Area Velocity vector
%W Characteristic waves vector - [forward; backward] waves
% R,lambda - right matrix, matrix of eigen vectors
%AFMh TVD edge Flux^hat
%WTMh edge characteristic wave ^hat
%r, phi - TVD parameter, flux corrector


N=ncell;
Q=[AREA; AU];
for i=2:N
    W(:,i)=R(:,:,i)\(Q(:,i)-Q(:,i-1));
end
W(:,1)=R(:,:,1)\(Q(:,1)-(3*Q(:,1)-Q(:,2))/2);
W(:,N+1)=R(:,:,N+1)\((3*Q(:,N)-Q(:,N-1))/2-Q(:,N-1));

r(1,1:N)=W(1,2:N+1)./W(1,1:N);
r(1,N+1)=((3*W(1,N+1)-W(1,N))/2)/W(1,N+1);

r(2,1:N)=W(2,2:N+1)./W(2,1:N);
r(2,N+1)=((3*W(2,N+1)-W(2,N))/2)/W(2,N+1);

if KTVD==1
    for i=1:N+1
        % Superbee
%         phi(1,i)=max([0, min(1,2*r(1,i)),min(2,r(1,i))]);
%         phi(2,i)=max([0, min(1,2*r(2,i)),min(2,r(2,i))]);
        
        % Van-Leer
        phi(1,i)=(r(1,i)+abs(r(1,i)))/(1+abs(r(1,i)));
        phi(2,i)=(r(2,i)+abs(r(2,i)))/(1+abs(r(2,i)));
        
        % Van Albada 1
%         phi(1,i)=(r(1,i)^2+r(1,i))/(1+r(1,i)^2);
%         phi(2,i)=(r(2,i)^2+r(2,i))/(1+r(2,i)^2);
        
        % Monotonized central
%         phi(1,i)=max(0,min( min(0.5*(1+r(1,i)),2),2*r(1,i)) );
%         phi(2,i)=max(0,min( min(0.5*(1+r(2,i)),2),2*r(2,i)) );
    end
else
    phi(1:2,1:N+1)=1;
end

for j=1:N+1
    WTMh(:,j)=phi(:,j).*W(:,j);
    AFMh(:,j)=0.5*R(:,:,j)*abs(Lambda(:,:,j))*(eye(2)-dtau/dx*abs(Lambda(:,:,j)))*WTMh(:,j);
end

end  %function TVDscalar


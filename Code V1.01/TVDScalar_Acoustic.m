function [AP,AU]=TVDScalar_Acoustic(AP,AU,cfl,rho, c, KTVD)
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


%Y-invariant
% N=length(Y);
cfl1=cfl*(1-cfl)/2;
Z=(AU-AP/rho/c);
Y=AU+AP/rho/c;
N=length(Y);
for i=2:N-1
    if i>2
        DYM=Y(i-1)-Y(i-2); %from page 2 theory look at r
        DY0=Y(i)  -Y(i-1);
        DYP=Y(i+1)-Y(i);
        %=====================
        if(DYM*DY0)<=0
            thm=0;
        else
            thm=DYM/DY0;
        end
        if(DYP*DY0)<=0
            thp=0;
        else
            thp=DY0/DYP;
        end
    else
        thp=0;
        thm=0;
    end
    
    %         fim=(thm+abs(thm))/(1+abs(thm));  %Van-Leer
    %         fip=(thm+abs(thp))/(1+abs(thp));
    %         fim=max(0,min( min(0.5*(1+thm),2),2*thm) );  %MC
    %         fip=max(0,min( min(0.5*(1+thp),2),2*thp) );  %MC
    fim=max([0, min(1,2*thm),min(2,thm)]) ;  %Superbee
    fip=max([0, min(1,2*thp),min(2,thp)]) ;  %Superbee
    if(KTVD==0)
        fip=1;fim=1;  %LW
    else
        if(KTVD==-1)
            fip=thp; fim=thm;  %B-W
        end
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% some type of flux here
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% some type of flux here    APT(2:N)=0.5*(AP(1:N-1)+AP(2:N))-0.5*dtau/h*(FLUXC1(2:N)-FLUXC1(1:N-1));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% some type of flux here    AUT(2:N)=0.5*(AU(1:N-1)+AU(2:N))-0.5*dtau/h*(FLUXC2(2:N)-FLUXC2(1:N-1));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% some type of flux here

    YN(i)=Y(i)-cfl*(Y(i)-Y(i-1))-cfl1*(fip*(Y(i+1)-Y(i))-fim*(Y(i)-Y(i-1)));

end
%     FLUX1=rhoc*c*AURIEM;
%     FLUX2=APRIEM/rho;
%     Y(1)= 
    Y(2:N-1)=YN(2:N-1);
%Z-invariant
%Z - invariant

for i=2:N-1
    if i>1 && i<N-2
    DZ0=Z(i)  -Z(i-1);
    DZP=Z(i+1)-Z(i);
    DZPP=Z(i+2)-Z(i+1);
    %=====================
    if(DZP*DZ0)<=0
        thm=0;
    else
        thm=DZP/DZ0;
    end
    if(DZP*DZPP)<=0
        thp=0;
    else
        thp=DZPP/DZP;
    end
       else
        thp=0;
        thm=0;
    end
    
    %         fim=(thm+abs(thm))/(1+abs(thm));  %Van-Leer
    %         fip=(thm+abs(thp))/(1+abs(thp));
    %         fim=max(0,min( min(0.5*(1+thm),2),2*thm) );  %MC
    %         fip=max(0,min( min(0.5*(1+thp),2),2*thp) );  %MC
    fim=max([0, min(1,2*thm),min(2,thm)]) ;  %Superbee
    fip=max([0, min(1,2*thp),min(2,thp)]) ;  %Superbee
    if(KTVD==0)
        fip=1;fim=1;  %LW
    else
        if(KTVD==-1)
            fip=thp; fim=thm;  %B-W
        end
    end
    ZN(i)=Z(i)+cfl*(Z(i+1)-Z(i))-cfl*(-cfl+1)/2*(fip*(Z(i+1)-Z(i))-fim*(Z(i)-Z(i-1)));
end %i
Z(2:N-1)=ZN(2:N-1);
%Z-invariant
%Z - invariant
    AU=0.5*(Y+Z);
    AP=0.5*rho*c.*(Y-Z);
    
end  %function TVDscalar


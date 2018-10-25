function [COMPLIANCEFINAL,RC,RP,RT]=RCRparameters(NVESSEL,LL,CMK,AREAZ,XX,NCELL,IVexit,NEXIT,RESIST,COMPLIANCE,Rho,Rmult,Cmult,correctcompliance,Rcorrect)
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

% for 26 we  have R total and R for each vessel CT is given
% for 55 and RT12 and Cj no need for calculation
% for 37 we  have RT12, Cj no need for calculation

%% Resistance Calculation
RTinv=0;
for jv=1:NEXIT
    RCd=0;
    RAd=0;

    % option1
        for NODE=1:NCELL(IVexit(jv))
            RAd=RAd+AREAZ(IVexit(jv),NODE)*XX(IVexit(jv),NODE);
            RCd=RCd+CMK(IVexit(jv),NODE)*XX(IVexit(jv),NODE);
        end
        RC(jv)=Rho*RCd/RAd;
    % option2
%         RC(jv)=Rho*CMK(IVexit(jv),NCELL(IVexit(jv)))/AREAZ(IVexit(jv),NCELL(IVexit(jv)));
    %
    if NVESSEL==26 && IVexit(jv)==20
        RC(jv)=Rcorrect*RC(jv);
    end
    RT(jv)=RESIST(jv)*Rmult;
    RP(jv)=RT(jv)-RC(jv);
    RTinv=RTinv+1/RT(jv);
end
Rtotal=1/RTinv;
% Rtotal=2.25*10^8; % total peripheral Resistances 37
% Rtotal=1e+8;  % total Peripheral Resistances 55


%% COMPLIANCES
if NVESSEL==26
% 26 vessel case only
    CC=0;
    Ad=0;
    Cd=0;
    for iv=1:NVESSEL
        for j=1:NCELL(iv)
            Ad=Ad+AREAZ(iv,j)*XX(iv,j);
            Cd=Cd+CMK(iv,j)*XX(iv,j);
        end
        Cseg=Ad/Rho/(Cd/LL(iv))^2; %m^3/Pa
        CC=CC+Cseg; %total conduit
    end
    CT=COMPLIANCE/Cmult; %1 number COMPLIANCE is CT
    CP=CT-CC;      %total peripheral
    COMPLIANCE=[CP*(69.1/100) CP*(19.3/100) CP*(5.2/100) CP*(6.4/100)];
    for kv=1:NEXIT
        COMPLIANCEFINAL(kv)=COMPLIANCE(kv)*(RT(kv))/RP(kv); %correction for 26 vessels case
    end  
else
% other cases
    for kv=1:NEXIT
        COMPLIANCEFINAL(kv)=COMPLIANCE(kv)/Cmult; %correction for 26 vessels case
    end
end


%% Correcting Compliance
% Correcting Compliance?
if NVESSEL~=26 && correctcompliance==1
    for kv=1:NEXIT
        COMPLIANCEFINAL(kv)=COMPLIANCEFINAL(kv)*(Rtotal)/RP(kv); %correction for 26 vessels case
    end
end

end
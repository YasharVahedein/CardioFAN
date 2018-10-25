function  [AP, AU]=RUN_NET_ACOU_INIT(SCALAR,ARRAY1D,AP,AU,AREAZ,CMK,XX,NODE_CONNECT,Pfit,Rmultiply,Cmultiply,Pout,acousticplot)
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
%APP,AUU - edge cell associated, APP(1:NVESSEL, 1:NCELLMAX+1)
%NODE_CONNECT - node connectivity matrix,
%NODE_CONNECT=[ALLNODES(1:NNODE),ivup,ivdn1,ivdn2], ivup, ivdn - global numbers
%ATIME(itau) - Array of time
%of upstream and doownstrean vessels

%Unpack SCALAR
%=============================================
rho=SCALAR(1);
KR=SCALAR(2);
NVESSEL=round(SCALAR(3)); NVESSEL=round(NVESSEL);
NNODE=round(SCALAR(4));
KTVD=SCALAR(5);
CFL=SCALAR(6);
dtau=SCALAR(7);
NTAU=SCALAR(8);
IVinlet=round(SCALAR(9));
NEXIT=round(SCALAR(10));
T=SCALAR(11);
TT=SCALAR(12);
TIME=SCALAR(13);
RCR=SCALAR(14);
Pzero=SCALAR(15);

%Unpack Arrays
%=============================================
n=1;m=0; m=m+NVESSEL;  NCELL= round(ARRAY1D(n:m));
n=m+1;   m=m+NVESSEL;  HMESH= ARRAY1D(n:m);
n=m+1;   m=m+NVESSEL;  LL=    ARRAY1D(n:m);
n=m+1;   m=m+NVESSEL;  ALF=   ARRAY1D(n:m);
n=m+1;   m=m+NEXIT;    IVEXIT= round(ARRAY1D(n:m));
n=m+1;   m=m+NEXIT;    RT=ARRAY1D(n:m);
n=m+1;   m=m+NEXIT;    RC=ARRAY1D(n:m);
n=m+1;   m=m+NEXIT;    RP=ARRAY1D(n:m);
n=m+1;   m=m+NEXIT;    COMPLIANCE=ARRAY1D(n:m);

%Define Max Cell number and initialize edge values of P and V
%=============================================
NCELLMAX=max(NCELL);
APP(1:NVESSEL,1:NCELLMAX)=AP(1:NVESSEL,1:NCELLMAX);
APP(1:NVESSEL,NCELLMAX+1)=AP(1:NVESSEL,NCELLMAX);
AUU(1:NVESSEL,1:NCELLMAX+1)=0;
AAA(1:NVESSEL,1:NCELLMAX)=AREAZ(1:NVESSEL,1:NCELLMAX);
AAA(1:NVESSEL,NCELLMAX+1)=AREAZ(1:NVESSEL,NCELLMAX);
AREA(1:NVESSEL,1:NCELLMAX)=AREAZ(1:NVESSEL,1:NCELLMAX);
%Initial P and V PLOT to SHOW INITIALIZATION in Acoustic code - Uncomment
%to see the initialization of accoustic
%=============================================
% for iv=1:NVESSEL
%     plot(XX(iv,:),AP(iv,:), 'linewidth',2);
%     grid on; title('INITIAL Pressure'); % ylim([1 6]);
%     figure
%     plot(XX(iv,1:NCELL(iv)),AU(iv,1:NCELL(iv)),'linewidth',2);
%     grid on; title('INITIAL Velocity');%ylim([-0.5 1.5]);%
% end

for i=1:NEXIT
    iv=IVEXIT(i);  nlast=NCELL(iv);
    AREABNp(iv)=AAA(iv,nlast);
    APBNp(iv)=AP(iv,nlast);
    AUBNp(iv)=AU(iv,nlast);
end

% Time loop Starts - 1/4 of the total time in Nonlinear Case (1 cycle) -
% Can be increased by changing NTAU/4
%=============================================
for itau=1:NTAU
    %==================
    
    t=dtau*itau;             %time flow
    ATIME(itau)=t;           %Array of time
    
    %Inlet Riemann's variables
    %==========================================================================================
    [ APP(IVinlet,1),  AUU(IVinlet,1) ] = INLET_Acoustic(AP(IVinlet,1), AU(IVinlet,1),AREA(IVinlet,1), rho*CMK(IVinlet,1), t-dtau/2,T,TT,Pfit);
    
    
    %==========================================================================================
    %Exit Riemann variables
    %======================================================================
    for i=1:NEXIT
        iv=IVEXIT(i);  nlast=NCELL(iv);  %exit vessel number and relating number of cells
        %         ref=(RESIST(i)-rho*CMK(iv)/(AREAZ(iv)))/(RESIST(i)+rho*CMK(iv)/(AREAZ(iv)));
        %==================================================================================
        % Complex
        [zU,AUBNp(iv), zP,APBNp(iv), zA,AREABNp(iv)]=... %added AREABNp, APNp, AUNp
            EXIT_Accoustic_Complex(AP(iv,nlast),APBNp(iv),AU(iv,nlast),AUBNp(iv),AREAZ(iv,nlast),...
            AREABNp(iv),rho,dtau,t,LL(iv),CMK(iv,nlast),RT(i),RC(i),RP(i),COMPLIANCE(i),Pout,RCR,Pzero);
        % Simple
        %         [zP,zU,zA]=EXIT_Acoustic(AP(iv,nlast),AU(iv,nlast),APP(iv,nlast+1),AUU(iv,nlast+1),AREAZ(iv,nlast),dtau,rho,CMK(iv,nlast),RESIST(i),COMPLIANCE(i));
        APP(iv,nlast+1)=zP;
        AUU(iv,nlast+1)=zU;
        AAA(iv,nlast+1)=zA;
        %==================================================================================
    end
    
    %Internal bifurcations Riemann variables
    %=====================================================================
    for inode=1:NNODE   %!!!make universal, avoiding inlet node and exit nodes
        %=========================================================================================
        ivup=NODE_CONNECT(inode,2);
        ivdn1=NODE_CONNECT(inode,3);
        ivdn2=NODE_CONNECT(inode,4);
        
        %Calculate P1N, U1N - pressure and velocity at the last cell of
        %parental vessel (number 1)
        %Calculate P21, U21 - pressure and velocity at the 1st cell of the
        %1st daughter vessel
        %Calculate P31, U31 - pressure and velocity at the 1st cell of the
        %2nd daughter vessel
        
        if(ivup~=0 && ivdn1~=0 )    %ivup=0 (inlet) and ivdn1=0 - exit - evaluated separately
            Nlast=NCELL(ivup);     %last cell at upstream vessel
            P1N=AP(ivup,Nlast);             U1N=AU(ivup,Nlast); %upstream centered values at the last cell
            roc1=rho*CMK(ivup,Nlast);       A1=AREA(ivup,Nlast);
            P21=AP(ivdn1,1);                U21=AU(ivdn1,1);    %downstream centered variable at the 1st cell
            roc2=rho*CMK(ivdn1,1);          A2=AREA(ivdn1,1);
            
            if(ivdn2~=0)
                P31=AP(ivdn2,1);              U31=AU(ivdn2,1);  %another      downstream local numeration 1st cell
                roc3=rho*CMK(ivdn2,1);        A3=AREA(ivdn2,1);
            else
                P31=0; U31=0; A3=0; roc3=10000;  %actually no bifurcation, byt Area change
            end
            
            [APP(ivup,Nlast+1),AUU(ivup,Nlast+1),APP(ivdn1,1),AUU(ivdn1,1),zP,zU]=...
                BIFURC_Acoustic(P1N,U1N,P21,U21,P31,U31,roc1,roc2,roc3,A1,A2,A3);
            if(ivdn2~=0)
                APP(ivdn2,1)=zP;
                AUU(ivdn2,1)=zU;
            end
        end
    end     %inode=1:NNODE
    
    %Riemann variables at internal cell facets for each vessel
    %=========================================================
    for iv=1:NVESSEL
        h=HMESH(iv);   ncell=NCELL(iv);
        
        for i=2:ncell   %i=1-the 1st cell &i=NCELL+1 are already evaluated
            roc=rho*CMK(iv,i);
            P1=AP(iv,i-1); P2=AP(iv,i); U1=AU(iv,i-1); U2=AU(iv,i);
            
            [APP(iv,i),AUU(iv,i)]=RIEMANN_Acoustic(P1,P2,U1,U2,roc);
            
        end
    end
    
    %Update at the next time moment cell centered values
    %============================================================
    for iv=1:NVESSEL
        h=HMESH(iv);  ncell=NCELL(iv);
        zU=dtau/(h*rho);  zP=dtau/h*(rho);
%         AU(iv,1)=AU(iv,1)-zU*diff(APP(iv,1:2));
%         AP(iv,1)=AP(iv,1)-zP*diff(AUU(iv,1:2).*CMK(iv,1:2).^2);
%         AU(iv,ncell)=AU(iv,ncell)-zU*diff(APP(iv,ncell:ncell+1));
%         AP(iv,ncell)=AP(iv,ncell)-zP*diff(AUU(iv,ncell:ncell+1).*CMK(iv,ncell1:ncell+1).^2);

        for i=1:ncell
            AU(iv,i)=AU(iv,i)-zU*(APP(iv,i+1)-APP(iv,i));
            AP(iv,i)=AP(iv,i)-zP*(AUU(iv,i+1)-AUU(iv,i))*CMK(iv,i)^2;  
        end
        
        [AP(iv,1:ncell),AU(iv,1:ncell)]=TVDScalar_Acoustic(AP(iv,1:ncell),AU(iv,1:ncell),CFL,rho,CMK(iv,1:ncell),KTVD);

        for j=1:ncell
            AREA(iv,j)=AREAZ(iv,j)*(((AP(iv,j)-Pzero)/(2*rho*(CMK(iv,j)^2)))+1)^2;
        end
        
        AU0(iv,itau)=AU(iv,1);
        AP0(iv,itau)=AP(iv,1);
        AA0(iv,itau)=AREA(iv,1);
        AUN(iv,itau)=AU(iv,round(ncell/2));
        APN(iv,itau)=AP(iv,round(ncell/2));
        AAN(iv,itau)=AREA(iv,round(ncell/2));

        %Animation of P or V as a function of space in vessel (vessel number)
        %============================================================
        %                         vesselnum=2;
        %                         if iv==vesselnum    %YS
        %                         figure (1)
        %                         plot(XX(iv,1:NCELL(iv)), AP(iv,1:NCELL(iv)),'LineWidth',3);
        %                         title('VESSEL-1 PRESSURE VS DISTANCE')
        %                         xlim([0 0.045])
        %                         ylim([-5000 12000])
        %                         ylabel('Pressure (Pa)') % y-axis label
        %                         figureHandle = gcf;
        %                         size = [200 200];
        %                         res = 300;
        %                         set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
        %                         set(findall(figureHandle,'type','text'),'fontSize',28, 'fontname', 'arial')
        %                         set(gcf, 'Color', 'w');
        %                         set(gca,'fontsize',28);
        %                         drawnow
        %                     end
        
    end
    
end

%Accoustic initialization in all vessels - uncomment for plots of
%initialization
%=========================================================
if acousticplot==1
    for iii=1:NVESSEL
        figure(iii+2)
        Casenum = int2str(iii);
        subplot(2,2,1)
        plot(ATIME, AU0(iii,:).*AA0(iii,:)*1000000,'Linewidth',2); title(['Acoustic Solution:Vessel ', Casenum, ' Inlet Velocity']);
        grid on
        
        subplot(2,2,2)
        plot(ATIME, ((AP0(iii,:))/1000),'Linewidth',2); title(['Acoustic Solution:Vessel ', Casenum, ' Inlet Pressure']);
        grid on
        
        subplot(2,2,3)
        plot(ATIME, AUN(iii,:).*AAN(iii,:)*1000000,'Linewidth',2); title(['Acoustic Solution:Vessel ', Casenum, ' Outlet Velocity']);
        grid on
        
        subplot(2,2,4)
        plot(ATIME, ((APN(iii,:))/1000),'Linewidth',2); title(['Acoustic Solution:Vessel ', Casenum, ' Outlet Pressure']);
        grid on
        hold on
    end
end
end




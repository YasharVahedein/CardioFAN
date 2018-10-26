function  RUN_NET_LW
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

%============NOT ACTIVE YET======================================
% parpool('local')      %Activate for parallel processing

%================================================================
% PARAMETER DESCRIPTION
%================================================================
%SCALAR  - array of all scalar inputs
%ARRAY1D - compound array of all 1D input arrays
%AP,AU   - 2D array - center cells 
%APB,AUB, AREAB - edge cell associated,
%FB,SB,HB - 
%NODE_CONNECT - node connectivity matrix,
%NODE_CONNECT=[ALLNODES(1:NNODE),ivup,ivdn1,ivdn2], ivup, ivdn - global numbers
% AT(itau) - Array of time

%%Clear 
close all; clc; 
%=============================================
%% READ INPUT FUNCTION
%=============================================
[SCALAR,ARRAY1D,PTTfunc,AP,AU,AREA,AREAZ,CMK,XX ,NODE_CONNECT,Pfit]=INPUT_TVD;

%=============================================
% Unpack SCALAR
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
a=SCALAR(14);
Pout=SCALAR(15);
RCR=SCALAR(16);
Pzero=SCALAR(17);

%=============================================
% Unpack Arrays
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

%=============================================
% Define Max Cell number
%=============================================
NCELLMAX=max(NCELL);

%=============================================%
% Initialize P, V and AREA in previous time-step for RCR TERMINAL boundary model
%=============================================%
for i=1:NEXIT
    iv=IVEXIT(i);  nlast=NCELL(iv);
    AREABNp(iv)=AREA(iv);
    APBNp(iv)=AP(iv,nlast);
    AUBNp(iv)=AU(iv,nlast);
end

%=============================================%
% Initialize first and last edges of each vessel + Fluxes and Variables
%=============================================%
for iv=1:NVESSEL 
    AUB1(iv)=AU(iv,1);              APB1(iv)=AP(iv,1);              AREAB1(iv)=AREAZ(iv,1);
    AUBN(iv)=AU(iv,NCELL(iv));      APBN(iv)=AP(iv,NCELL(iv));      AREABN(iv)=AREAZ(iv,NCELL(iv)); 
end
    Q2(1:2,1:NCELLMAX,1:NVESSEL)=0;
    AF2(1:2,1:NCELLMAX,1:NVESSEL)=0;
    S2(1:2,1:NCELLMAX,1:NVESSEL)=0;
    H(1:2,1:2,1:NCELLMAX,1:NVESSEL)=0;
    B(1:2,1:2,1:NCELLMAX,1:NVESSEL)=0;
    
    
  
    
    
    
%=============================================
%% Time loop Starts
%=============================================
AT(1)=0;  t=0;          %Initialize time and its array

for itau=1:100000000    %Run Until reahcing TIME
    
    
    %=============================================
    %% INLET & OUTLET
    %=============================================
    %Inlet Edge variables
    [ AUB1(IVinlet),  APB1(IVinlet), AREAB1(IVinlet), FB1(:,IVinlet), HB1(:,:,IVinlet), SB1(:,IVinlet)] =...
        INLET_HYPER(AU(IVinlet,1), AREA(IVinlet,1),AREAZ(IVinlet,1), rho,dtau, t+dtau/2,TT,CMK(IVinlet,1),KR,Pfit,a,Pzero);

    %Exit Edge variables
    for i=1:NEXIT
        iv=IVEXIT(i);  nlast=NCELL(iv);  %exit vessel number and its number of cells
        [AUBN(iv),AUBNp(iv), APBN(iv),APBNp(iv), AREABN(iv),AREABNp(iv), FBN(:,iv), HBN(:,:,iv), SBN(:,iv)]=... %added AREABNp, APNp, AUNp
            EXIT_HYPER(AP(iv,nlast),APBNp(iv),AU(iv,nlast),AUBNp(iv),AREAZ(iv,nlast),AREA(iv,nlast),AREABNp(iv)...
            ,rho,dtau,t,LL(iv),CMK(iv,nlast),KR,RT(i),RC(i),RP(i),COMPLIANCE(i),a,Pout,RCR,Pzero);
    end

    %=============================================
    % Time advancing
    %=============================================
    t=t+dtau;                                           % time
    AT(itau)=t;                                         %ARRAY of Time
    if(t>TIME) || imag(t)~=0
        AT=AT(1:itau-1);                                %Accumulate time array and break when it ends
        break
    end
    
    
    %% JOINT BOUNDARY CONDITIONS
    %=====================================================================
    for inode=1:NNODE   % Starting from Node 1 to NNODE
    %=====================================================================
        ivup=NODE_CONNECT(inode,2);    %Upstream Vessel 
        ivdn1=NODE_CONNECT(inode,3);   %Downstream Vessel1
        ivdn2=NODE_CONNECT(inode,4);   %Downstream Vessel2
        
        %Calculate P1N, U1N - pressure and velocity at the last cell of
        %parental vessel (number 1)
        %Calculate P21, U21 - pressure and velocity at the 1st cell of the
        %1st daughter vessel
        %Calculate P31, U31 - pressure and velocity at the 1st cell of the
        %2nd daughter vessel
        
        if(ivup~=0 && ivdn1~=0 )    %ivup=0 (inlet) and ivdn1=0 - exit - evaluated separately
            Nlast=NCELL(ivup);     %last cell at upstream vessel
            P1N=AP(ivup,Nlast);  U1N=AU(ivup,Nlast); %upstream centered values at the last cell
            roc1=rho*CMK(ivup,Nlast);     A1=AREA(ivup,Nlast);   A01N=AREAZ(ivup,Nlast);   cmk1N=CMK(ivup,Nlast);  kr1=KR;
            P21=AP(ivdn1,1);     U21=AU(ivdn1,1);    %downstream centered variable at the 1st cell
            roc2=rho*CMK(ivdn1,1);    A2=AREA(ivdn1,1);      A021=AREAZ(ivdn1,1);       cmk21=CMK(ivdn1,1);  kr2=KR;
            
            if(ivdn2~=0)
                P31=AP(ivdn2,1);     U31=AU(ivdn2,1);  %another      downstream local numeration 1st cell
                roc3=rho*CMK(ivdn2,1);    A3=AREA(ivdn2,1);      A031=AREAZ(ivdn2,1);    cmk31=CMK(ivdn2,1); kr3=KR;
            else
                P31=0; U31=0;  roc3=10000; A3=0; A031=0; cmk31=0; kr3=KR; %YS - maybe AREAZ(ivdn2); %actually no bifurcation, byt Area change
            end
            
            %Calculates the last and first edge centered cell values of parent and daughter, respectively
            %=====================================================================
            [APBN(ivup),AUBN(ivup), AREABN(ivup), APB1(ivdn1),AUB1(ivdn1), AREAB1(ivdn1),zP,zU,zA, ...
                FBN(:,ivup), FB1(:,ivdn1),zFB, HBN(:,:,ivup),HB1(:,:,ivdn1),zHB,SBN(:,ivup),SB1(:,ivdn1),zSB]=...
                BIFURC_HYPER(P1N,U1N,P21,U21,P31,U31,roc1,roc2,roc3,A1,A01N,A2,A021,A3,A031,cmk1N,cmk21,cmk31,kr1,kr2,kr3,rho,a,Pzero);
            if(ivdn2~=0)
                APB1(ivdn2)=zP;
                AUB1(ivdn2)=zU;
                AREAB1(ivdn2)=zA;
                FB1(:,ivdn2)=zFB;
                HB1(:,:,ivdn2)=zHB;
                SB1(:,ivdn2)=zSB;
            end
        end
    end     %inode=1:NNODE
    
     %Update at the next time moment cell centered values
    %============================================================
    for iv=1:NVESSEL
        dx=HMESH(iv);   ncell=NCELL(iv); 
        %=============================================
        %Cell centered vectors and matrices
        %====================================
        for k=1:ncell
            Q2(:,k,iv)=[AREA(iv,k); AU(iv,k)];
            AF2(:,k,iv)=[AU(iv,k)*AREA(iv,k); AU(iv,k)^2/2+AP(iv,k)/rho];
            S2(:,k,iv)=[0; KR*AU(iv,k)/AREA(iv,k)];
            H(:,:,k,iv)=[AU(iv,k) AREA(iv,k); (CMK(iv,k)^2)/sqrt(AREAZ(iv,k)*AREA(iv,k))*exp(a*(sqrt(AREA(iv,k)/AREAZ(iv,k))-1)^2)*(1+2*a*(sqrt(AREA(iv,k)/AREAZ(iv,k))-1)^2) AU(iv,k)]; %--VAR AREA--%
            B(:,:,k,iv)=KR*[0 0; -AU(iv,k)/AREA(iv,k)^2 1/AREA(iv,k)]; 
%             AQ2(:,k,iv,itau)=Q2(:,k,iv);
        end
        %Averaging for cell edges -                                         %--VAR AREA--%  define S ave 
        %===============================================
        HAV(:,:,1,iv)  =[AUB1(iv) AREAB1(iv); (CMK(iv,1)^2)/sqrt(AREAZ(iv,1)*AREAB1(iv))*exp(a*(sqrt(AREAB1(iv)/AREAZ(iv,1))-1)^2)*(1+2*a*(sqrt(AREAB1(iv)/AREAZ(iv,1))-1)^2) AUB1(iv)];
        HAV(:,:,ncell+1,iv)=[AUBN(iv) AREABN(iv); (CMK(iv,ncell)^2)/sqrt(AREAZ(iv,ncell)*AREABN(iv))*exp(a*(sqrt(AREABN(iv)/AREAZ(iv,ncell))-1)^2)*(1+2*a*(sqrt(AREABN(iv)/AREAZ(iv,ncell))-1)^2) AUBN(iv)]; %--VAR AREA--%  think about AREAZ
        for k=2:ncell
            HAV(:,:,k,iv)=0.5*(H(:,:,k,iv)+H(:,:,k-1,iv));
        end
        %====================================================
        %Update centered cell values                                         % check s values %--VAR AREA--%
        %===================================================
        for k=1:ncell
            if(k==1)% && iv==1)
                H0=2*HB1(:,:,iv)-H(:,:,1,iv);
                S20=2*SB1(:,iv)-S2(:,1,iv);
                AF20=2*FB1(:,iv)-AF2(:,1,iv); 
            else
                H0=H(:,:,k-1,iv);
                S20=S2(:,k-1,iv);
                AF20=AF2(:,k-1,iv);
            end
            if(k==ncell)% && iv>1)
                HNp1=2*HBN(:,:,iv)-H(:,:,k,iv); 
                S2N=2*SBN(:,iv)-S2(:,k,iv);     
                AF2N=2*FBN(:,iv)-AF2(:,k,iv);   
            else
                HNp1=H(:,:,k+1,iv);
                S2N=S2(:,k+1,iv);
                AF2N=AF2(:,k+1,iv);
            end
            
            DF=AF2N-AF20; % Fi+1-Fi-1
            DFP=AF2N-AF2(:,k,iv); % Fi+1-Fi
            DFM=AF2(:,k,iv)-AF20; % Fi-Fi-1
            
            DQ=(eye(2)+dtau/2*B(:,:,k,iv))*(S2(:,k,iv)-DF/(2*dx));
            D2Q=((HAV(:,:,k+1,iv)*(DFP)-HAV(:,:,k,iv)*(DFM))/dx^2-((HNp1*S2N-H0*S20)/(2*dx)));

            Q2(:,k,iv)=Q2(:,k,iv)+dtau*DQ+(dtau^2)/2*D2Q;
            AREA(iv,k)=Q2(1,k,iv);
            AU(iv,k)=Q2(2,k,iv);
            AP(iv,k)=2*rho*CMK(iv,k)^2*(sqrt(AREA(iv,k)/AREAZ(iv,k))-1)*exp(a*(sqrt(AREA(iv,k)/AREAZ(iv,k))-1)^2)+Pzero; %% new change in AP calculation
        end  %k loop
        
        %============================================================
        %% Animation of P or V as a function of space in desired vessel (vessel number)
        %============================================================
        %           vesselnum=1;
        %           if iv==vesselnum    %YS
        %               figure (1)
        %               AUdisp=Q2(2,1:ncell,iv);
        %               plot(XX(1,1:ncell), AUdisp(1,1:ncell),'LineWidth',3);
        %               title('VESSEL-1 Velo VS DISTANCE')
        % %               xlim([0 0.04])
        % %               ylim([0 10])
        %               ylabel('Velo (mm/s)') % y-axis label
        %               figureHandle = gcf;
        %               size = [200 200];
        %               res = 300;
        %               set(gcf,'paperunits','inches','paperposition',[0 0 size/res]);
        %               set(findall(figureHandle,'type','text'),'fontSize',28, 'fontname', 'arial')
        %               set(gcf, 'Color', 'w');
        %               set(gca,'fontsize',28);
        %               drawnow
        %           end
        
        %===============================
        %% Collect time dependent arrays for 1st and middle cell at each vessel
        %===============================
        AREA0(iv,itau)=AREA(iv,1);
        AU0(iv,itau)=AU(iv,1);
        AP0(iv,itau)=AP(iv,1);
        AREAN(iv,itau)=AREA(iv,round(ncell/2));
        AUN(iv,itau)=AU(iv,round(ncell/2));
        APN(iv,itau)=AP(iv,round(ncell/2));
        if PTTfunc(1)==1 && iv<=PTTfunc(5)
            PWV2(iv,1:ncell,itau)=AU(iv,1:ncell)+CMK(iv,1:ncell).*(exp(a*((sqrt(AREA(iv,1:ncell)./AREAZ(iv,1:ncell))-1).^2)/2)).*sqrt(sqrt(AREA(iv,1:ncell)./AREAZ(iv,1:ncell)).*(1+2*a*(sqrt(AREA(iv,1:ncell)./AREAZ(iv,1:ncell))-1).^2));
        end
    end
    
    %====================================================
    %% CFL based time-step size calculation
    %====================================================
    for iv=1:NVESSEL
        vmax=max(AU(iv,:)+CMK(iv,:).*sqrt(sqrt(AREA(iv,:)./AREAZ(iv,:)).*exp(a*(sqrt(AREA(iv,:)./AREAZ(iv,:))-1).^2).*(1+2*a*(sqrt(AREA(iv,:)./AREAZ(iv,:))-1).^2)));            % Find maximum eigenvalue in center variables in all vessels
        vmax=max(vmax,AUB1(iv)+CMK(iv,1).*sqrt(sqrt(AREAB1(iv)/AREAZ(iv,1))*exp(a*(sqrt(AREAB1(iv)/AREAZ(iv,1))-1)^2)*(1+2*a*(sqrt(AREAB1(iv)/AREAZ(iv,1))-1)^2))); % Find maximum eigenvalue in edge variables or with UFUN(Umax,t,T,TT)+... in all vessels
    end
    vmax=max(vmax(1,:));                                % Maximum Velocity in th whole Network
    dtau=(CFL*min(HMESH)/vmax);                         % Find dtau
end  %END time loop -itau loop





%===============================
%% PTT CALCULATION
%===============================
if PTTfunc(1)==1
    PTTcalc(PWV2,AT,HMESH,NCELL,PTTfunc,AU0,AREA0);
end


%============================================================
%% Final plot of P and V in first and mid cell at each vessel
%============================================================
for iii=1:NVESSEL
%         if iii==2 || iii==10 || iii==13 || iii==17 || iii==19 || iii==24
% %         if iii==19 || iii==24 
        figure(iii+NVESSEL+2)
        Casenum = int2str(iii);
        subplot(2,2,1)
        plot(AT, AU0(iii,:).*AREA0(iii,:)*1000000,'k','Linewidth',1); title(['Vessel ', Casenum, ' Inlet Flow Rate']); % was 7.4506e-04
        grid off
        hold on
        subplot(2,2,2)
        plot(AT, ((AP0(iii,:))/1000),'k','Linewidth',1); title(['Vessel ', Casenum, ' Inlet Pressure (kPa)']);
        grid off 
        hold on
        subplot(2,2,3)
        plot(AT, AUN(iii,:).*AREAN(iii,:)*1000000,'k','Linewidth',1); title(['Vessel ', Casenum, ' Mid-Seg Flow Rate']); % was 7.4506e-04
        grid off 
        hold on
        subplot(2,2,4)
        plot(AT, ((APN(iii,:))/1000),'k','Linewidth',1); title(['Vessel ', Casenum, ' Mid-Seg Pressure (kPa)']);
        grid off 
        hold on
        
        %calculate the vessel area
        figure(iii+NVESSEL+2+26)
        Casenum = int2str(iii);
        subplot(2,1,1)
        plot(AT, AREA0(iii,:)*10000,'k','Linewidth',1); title(['Vessel ', Casenum, ' Inlet AREA (cm^2)']); % was 7.4506e-04
        grid off
        hold on
        subplot(2,1,2)
        plot(AT, AREAN(iii,:)*10000,'k','Linewidth',1); title(['Vessel ', Casenum, ' Mid-Seg AREA (cm^2)']); % was 7.4506e-04
        grid off 
        hold on
%         end
%         fig=figure(iii+NVESSEL);
%         saveas(fig,sprintf('FIG%d.fig',iii))
end
end

    %===============================
    %Final stacked plots of P and V in first and mid cell of select vessels
    %change vessel numbers for your arbitrary plots
    %===============================
%     figure(4)
%     subplot(2,2,1)
%     plot(AT, AU0(1,:).*AREA0(1,:)*1000000,'Linewidth',2); title(['Nonlinear Solution:Vessel ', '1', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(1,:))/133.3)),'Linewidth',2); title(['Nonlinear Solution:Vessel ', '1', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(1,:).*AREAN(1,:)*1000000,'Linewidth',2); title(['Nonlinear Solution:Vessel ', '1', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(1,:))/133.3)),'Linewidth',2); title(['Nonlinear Solution:Vessel ', '1', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(14,:).*AREA0(14,:)*1000000,'Linewidth',4); title(['Nonlinear Solution:Vessel ', '14', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(14,:))/133.3)),'Linewidth',4); title(['Nonlinear Solution:Vessel ', '14', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(14,:).*AREAN(14,:)*1000000,'Linewidth',4); title(['Nonlinear Solution:Vessel ', '14', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(14,:))/133.3)),'Linewidth',4); title(['Nonlinear Solution:Vessel ', '14', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(27,:).*AREA0(27,:)*1000000,'-'); title(['Nonlinear Solution:Vessel ', '27', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(27,:))/133.3)),'-'); title(['Nonlinear Solution:Vessel ', '27', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(27,:).*AREAN(27,:)*1000000,'-'); title(['Nonlinear Solution:Vessel ', '27', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(27,:))/133.3)),'-'); title(['Nonlinear Solution:Vessel ', '27', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(39,:).*AREA0(39,:)*1000000,'--'); title(['Nonlinear Solution:Vessel ', '39', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(39,:))/133.3)),'--'); title(['Nonlinear Solution:Vessel ', '39', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,3)
%     plot(AT, AUN(39,:).*AREAN(39,:)*1000000,'--'); title(['Nonlinear Solution:Vessel ', '39', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(39,:))/133.3)),'--'); title(['Nonlinear Solution:Vessel ', '39', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     
%     
%     figure(5)
%     subplot(2,2,1)
%     plot(AT, AU0(15,:).*AREA0(15,:)*1000000,'Linewidth',2); title(['Nonlinear Solution:Vessel ', '15', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(15,:))/133.3)),'Linewidth',2); title(['Nonlinear Solution:Vessel ', '15', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(15,:).*AREAN(15,:)*1000000,'Linewidth',2); title(['Nonlinear Solution:Vessel ', '15', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(15,:))/133.3)),'Linewidth',2); title(['Nonlinear Solution:Vessel ', '15', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(21,:).*AREA0(21,:)*1000000,'Linewidth',4); title(['Nonlinear Solution:Vessel ', '21', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(21,:))/133.3)),'Linewidth',4); title(['Nonlinear Solution:Vessel ', '21', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(21,:).*AREAN(21,:)*1000000,'Linewidth',4); title(['Nonlinear Solution:Vessel ', '21', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(21,:))/133.3)),'Linewidth',4); title(['Nonlinear Solution:Vessel ', '21', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(38,:).*AREA0(38,:)*1000000,'-'); title(['Nonlinear Solution:Vessel ', '38', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(38,:))/133.3)),'-'); title(['Nonlinear Solution:Vessel ', '38', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,3)
%     plot(AT, AUN(38,:).*AREAN(38,:)*1000000,'-'); title(['Nonlinear Solution:Vessel ', '38', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(38,:))/133.3)),'-'); title(['Nonlinear Solution:Vessel ', '38', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold on
%     
%     
%     subplot(2,2,1)
%     plot(AT, AU0(46,:).*AREA0(46,:)*1000000,'--'); title(['Nonlinear Solution:Vessel ', '46', ' Inlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,2)
%     plot(AT, (((AP0(46,:))/133.3)),'--'); title(['Nonlinear Solution:Vessel ', '46', ' Inlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,3)
%     plot(AT, AUN(46,:).*AREAN(46,:)*1000000,'--'); title(['Nonlinear Solution:Vessel ', '46', ' Outlet Velocity']); xlabel('Time (7th cycle, 1s Interval)'); ylabel('Flow Rate (mL/s)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
%     subplot(2,2,4)
%     plot(AT, (((APN(46,:))/133.3)),'--'); title(['Nonlinear Solution:Vessel ', '46', ' Outlet Pressure']);      xlabel('Time (7th cycle, 1s Interval)'); ylabel('Pressure (mmHg)');
%     xlim([7 8]);
%     grid on
%     hold off
%     
% end

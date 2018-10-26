function PTTcalc(PWV2,AT,HMESH,NCELL,PTTfunc,AU0,AREA0)

STARTtime=PTTfunc(2);
FINISHtime=PTTfunc(3);
pttSTAvessel1=PTTfunc(4);
pttFINvessel2=PTTfunc(5);
LLLmax=PTTfunc(6);
% STARTtime=1;
% FINISHtime=5000;
PWVVVmin=50;
% PWV2(1,:,itau)=AU0(:,itau)+CMK(1,:)*sqrt(sqrt(AREA0(:,itau)/AREAZ(1,:)));


for Stime=STARTtime:FINISHtime
    
    for iv=pttSTAvessel1:pttFINvessel2
        for cell=1:NCELL(iv)
            PWVmin=min(min(PWV2(iv,cell,1+Stime:Stime+5000))); % was 1,:,:
            PWVVVmin=min(PWVmin,PWVVVmin);
        end
    end
    
    TTTmax=LLLmax/PWVVVmin;
    % TTTmax=sum(Tmax);
    if Stime==0
        AT1(1)=0;
        AT1(2:length(AT)+1)=AT;
    else
        AT1(1)=AT(Stime);
        AT1(2:Stime+5000-Stime+1)=AT(1+Stime:Stime+5000);
    end
    
    for ii=1:length(AT1)
        if AT1(ii)>=TTTmax+AT1(1)
            imax=ii;
            break
        end
    end
    
    dx=HMESH(pttSTAvessel1);
    %Accumulate PTT by subsegments
    PTT=dx/PWV2(pttSTAvessel1,1,1+Stime);  %time required to pass cell 1
    for cellnum=2:NCELL(pttSTAvessel1)
        %basic PWV array in time at i-th cell
        PWVbas(1:imax)=PWV2(2,cellnum,1+Stime:imax+Stime);
        PWVi = interp1(AT1(1:imax),PWVbas(1:imax),PTT+AT1(1),'linear'); %PWV at the center of next cell
        %time required to pass current cell
        PTT=PTT+dx/PWVi;
    end
    
    for iv=pttSTAvessel1+1:pttFINvessel2
        for cellnum=1:NCELL(iv)
            dx=HMESH(iv);
            %basic PWV array in time at i-th cell
            PWVbas(1:imax)=PWV2(iv,cellnum,1+Stime:imax+Stime);
            PWVi = interp1(AT1(1:imax),PWVbas(1:imax),PTT+AT1(1),'linear'); %PWV at the center of next cell
            %time required to pass current cell
            PTT=PTT+dx/PWVi;
        end
    end
    
    PTT;
    APTT(Stime)=PTT;
    %iQ loop
end

figure(5)
casenum1=int2str(pttSTAvessel1);
casenum2=int2str(pttFINvessel2);
hold on
subplot(2,1,1)
plot(AT(STARTtime:FINISHtime), APTT(STARTtime:FINISHtime),'k','Linewidth',3); title(['Vessel ',casenum1,' to ',casenum2, ' PTT as a Function of starting point']); % was 7.4506e-04
xlabel('Start Time of PTT Calculation'); ylabel('Pulse Transit Time (s)');
dim = [.2 .5 .3 .3];
AVE=num2str(mean(APTT(STARTtime:FINISHtime)));
STD=num2str(std(APTT(STARTtime:FINISHtime)));
MIN=num2str(min(APTT(STARTtime:FINISHtime)));
MAX=num2str(max(APTT(STARTtime:FINISHtime)));
str = {'PTT','Mean=',AVE,'STDev=',STD,'Min=',MIN,'Max=',MAX};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
grid off
hold on
subplot(2,1,2)
plot(AT(STARTtime:FINISHtime), AU0(pttSTAvessel1,STARTtime:FINISHtime).*AREA0(pttSTAvessel1,STARTtime:FINISHtime)*1000000,'k','Linewidth',1); title(['Vessel ', casenum1, ' Inlet Flow Rate']); % was 7.4506e-04
grid off
hold on
        % figure(5)
% plot(2*2000*(count/20)*0.25/pi,APTT,'gs','LineWidth',10);
% xlabel('Stroke Volume (ml/s)'); ylabel('Pulse Transit Time (s)')
% title('PTT as a Function of Stroke Volume')
% hold on
% grid on

end
function [A,AU1,Pfit]=FLOWINPUT(TIME,Period,A0,N)
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

if N==26
A=[0	-1.62455	;
0.002   -3   ;
0.005   -2   ;
0.008   4   ;
0.010   10  ;
0.018   18  ;
0.026	26.1228	;
0.052	166.67049	;
0.078	381.79218	;
0.104	460.17181	;
0.13	476.38968	;
0.156	434.53793	;
0.18200002	386.89871	;
0.208	334.58829	;
0.234	273.63593	;
0.26	214.64903	;
0.286	167.85107	;
0.312	111.42528	;
0.338	17.46402	;
0.36400003	-28.61443	;
0.39000003	-8.1915	;
0.416	3.00541	;
0.442	18.71727	;
0.468	33.14891	;
0.49400003	39.66362	;
0.52	36.8949	;
0.546	23.39226	;
0.572	16.76637	;
0.598	10.61592	;
0.624	9.46768	;
0.65000006	9.95389	;
0.676	8.48257	;
0.702	10.17103	;
0.72800006	8.56073	;
0.754	10.15691	;
0.78000006	11.3166	;
0.806	9.61741	;
0.832	5.31422	;
0.85800006	5.21226	;
0.884	5.12412	;
0.91	7.42971	;
0.936	9.07686	;
0.962	10.64409	;
0.979	2.26921	;
1.00    -1;
1.01400006	+2.67516];


Pfit= fit( A(1:46,1),A(1:46,2), 'smoothingspline');

elseif N==55 
    
%% 55
A=[0        6.55058;
0.01221     39.3864;
0.019536	88.6401;
0.028083	137.894;
0.035409	200.829;
0.040293	256.468;
0.047619	307.546;
0.0549451	353.151;
0.0610501	389.635;
0.0695971	425.207;
0.0830281	458.955;
0.0976801	486.318;
0.112332	500.912;
0.125763	507.297;
0.147741	495.439;
0.163614	472.637;
0.178266	438.889;
0.192918	413.35;
0.208791	379.602;
0.21978     337.645;
0.23199     302.073;
0.241758	259.204;
0.253968	221.808;
0.26862     188.972;
0.28083     154.312;
0.294261	104.146;
0.30525     38.4743;
0.310134	-4.39469;
0.315018	-30.8458;
0.321123	-43.6153;
0.330891	-14.4279;
0.338217	-4.39469;
0.349206	3.81426;
0.371184	7.46269;
0.401709	5.63847;
0.579976	8.37479;
0.65812     7.46269;
0.775336	9.2869;
0.885226	1.07794;
0.95        8.37479];

% Pfit=polyfit(A(1:56,1),A(1:56,2),17);
Pfit= fit( A(1:40,1),A(1:40,2), 'smoothingspline');


elseif N==37
A=[0.00E+00	9.259259
0.01024511	24.080004
0.014723141	39.823288
0.016944997	59.269
0.02141136	82.88265
0.02588184	103.71853
0.029606897	121.31324
0.034095224	130.11208
0.043839954	135.67314
0.04907439	147.25018
0.0550467	167.16098
0.063254595	197.25822
0.07147965	215.7814
0.08197323	222.26881
0.09921543	231.07486
0.10670948	236.63464
0.12471289	238.0337
0.14345418	247.76651
0.1592179	241.75688
0.17572017	243.61806
0.19525973	220.94388
0.2095271	212.15564
0.22755043	200.12877
0.26442966	120.51989
0.31562296	0.6412863
0.32540062	-16.019875
0.34041274	-21.103994
0.35019523	-41.005898
0.36453193	-96.553406
0.3728064	-111.363556
0.38629746	-102.55964
0.39376405	-78.48133
0.4049605	-40.04905
0.41469082	-24.765766
0.42520294	-30.778355
0.43646663	-37.716446
0.44772282	-39.561943
0.45821846	-34.46342
0.470973	-34.91918
0.48448056	-37.22637
0.51523924	-36.746044
0.5489858	-27.467724
0.5722339	-21.436075
0.6172421	-17.706957
0.63824713	-16.76917
0.66524506	-9.809478
0.693753	-9.33042
0.71774864	-1.4464946
0.7334822	12.914252
0.74545765	31.902512
0.75894874	40.706432
0.77246106	35.1585
0.7807307	23.589085
0.7950043	10.634171
0.807765	6.011741
0.82	9.259259];

% Pfit=polyfit(A(1:56,1),A(1:56,2),17);
Pfit= fit( A(1:56,1),A(1:56,2), 'smoothingspline');
end



%=========================================================
%% plot Inlet Velocity (UFUN) as a function of time
%=========================================================
At=linspace(0,TIME,500);                      %Time serries
tmax=3.0;                                     %to plot within tmax
t=0;                                          %current time

for i=1:2000
    At(i)=t;
    AU1(i)=UFUN_NONLIN(Pfit,At(i),Period,A0); % Inlet Velocity Function
    t=t+tmax/2000;                            % current time update
end
figure(1)
plot(At,AU1.*A0*1000000,'linewidth',2);xlabel('Time (sec)'); ylabel('U');
title('Inlet Velocity Time History'); grid on
end

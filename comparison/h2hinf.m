clc
clear 
close all

load comparison

AL=A;
BL=B;
BLw=H;
CL=C;
D=[0 0];

Bw=H;

%% Wadditive
syms r
dum1=expand((r+500)^7);
dum2=expand((r+700)^7);
[c1,t1] = coeffs(dum1);
[c2,t2] = coeffs(dum2);

l1=double(c1);
l2=double(c2);

WR=tf(l1,l2)*eye(2);

%% Wactuator
dum5=expand((r+1500)^1);
dum6=expand((r+50)^1);

[c5,t5] = coeffs(dum5);
[c6,t6] = coeffs(dum6);

l5=double(c5);
l6=double(c6);

Wact=tf(l5,l6)*eye(2);
% Wact=eye(2);


%% WPerformance
dum3=expand((r+1500)^1);
dum4=expand((r+50)^1);

[c3,t3] = coeffs(dum3);
[c4,t4] = coeffs(dum4);

l3=double(c3);
l4=double(c4);

Wperf=tf(l3,l4);
% Wperf=1;


%% Wdist
dum9=expand((r+200)^1);
dum10=expand((r+1)^1);

[c9,t9] = coeffs(dum9);
[c10,t10] = coeffs(dum10);

l9=double(c9);
l10=double(c10);

% Wdist=tf(l9,l10);
% Wdist=tf(l5,l6);
Wdist=1;

%%
% S1=ss(AL,BLw,CL,0);
% S2=ss(A,Bw,C,0);
% 
% % figure;
% % bodemag(S1,S2)
% 
% % (P_ho - P_Es)/P_Es < W(s)
% Upper_Limit=S1*(1+WR);
% Lower_Limit=S1;
% 
% figure
% bodemag(S1,'b',Upper_Limit,'r--',Lower_Limit,'k')
% 
% figure
% bodemag(S2,'k',WR,'r--')
% 
% figure
% bodemag(S2,'k',Wperf,'r--')
% 
% figure
% bodemag(S2,'k',Wact,'r--')
% 
% figure
% bodemag(S2,'k',Wdist,'r--')

%%
% AL=A;BL=B;BLw=Bw;CL=C;
sys1=ss(AL,BLw,CL,0);

%%%%%% Augmented Plant

sys2=ss(AL,[BLw BL],CL,[0 0 0]);

systemnames = 'sys2 WR Wperf Wact Wdist'; 
inputvar = '[dist;fs1;fs2]'; 
outputvar = '[WR;Wperf;Wact;Wperf]'; 
input_to_sys2 = '[Wdist;fs1;fs2]'; 
input_to_WR = '[fs1;fs2]'; 
input_to_Wperf = '[sys2]'; 
input_to_Wact = '[fs1;fs2]';
input_to_Wdist = '[dist]';
sys3 = sysic; 
%% Hinf
[K3,~,GAM] = hinfsyn(sys3,1,2);

%% H2Hinf
%%%% H2/Hinf Controller Design
r=[2 1 2];
obj=[10000000 0 0 0];


tt = pi/3;
bb = -1000;
region=[cos(tt) sin(tt) -2*bb*cos(tt) 0;-sin(tt) cos(tt) 0 -2*bb*cos(tt)];
% region=[];

[gopt,h2opt,K,R,S] = hinfmix(lti2mat(sys3),r,obj,region);

Kmix=mat2lti(K);

AC=Kmix.A;
BC=Kmix.B;
CC=Kmix.C;
DC=Kmix.D;
save('AC.mat','AC');
save('BC.mat','BC');
save('CC.mat','CC');
save('DC.mat','DC');
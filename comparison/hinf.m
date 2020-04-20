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
% bodemag(S2,'k',Wdist,'r--')

%%
% AL=A;BL=B;BLw=Bw;CL=C;
sys1=ss(AL,BLw,CL,0);

%%%%%% Augmented Plant

sys2=ss(AL,[BLw BL],CL,[0 0 0]);

systemnames = 'sys2 WR Wdist'; 
inputvar = '[dist;fs1;fs2]'; 
outputvar = '[sys2]'; 
input_to_sys2 = '[Wdist;fs1;fs2]'; 
input_to_WR = '[fs1;fs2]'; 
input_to_Wdist = '[dist]';
sys3 = sysic; 
%% Hinf
[K3,~,GAM] = hinfsyn(sys3,1,2);

AC=K3.A;
BC=K3.B;
CC=K3.C;
DC=K3.D;
save('AC.mat','AC');
save('BC.mat','BC');
save('CC.mat','CC');
save('DC.mat','DC');
clc
clear

A=[-209.443458108676,450.728478793718,0,0,0,0;-450.728478793718,-209.443458108676,0,0,0,0;0,0,-2.25680753735172,467.128469627320,0,0;0,0,-467.128469627320,-2.25680753735172,0,0;0,0,0,0,-0.460680331454473,87.5968670596031;0,0,0,0,-87.5968670596031,-0.460680331454473];
B=[-9.64474367821996,2.65660382188316;-7.76311006325836,-5.61613187293423;-0.117927290005266,0.0527206875641465;0.215603654718808,-0.128623447556356;0.0737096142170821,0.0391127484963023;-0.0970049785583168,-0.0786358112919118];
C=[1.88477251668799,-1.14284793707382,5.77230310808448,47.9525873535994,-5.11774457916239,-49.0010462763406];
H=-B(:,1);

dum1=size(A);
dum2=size(B);
dum3=size(H);
dum4=size(C);
n=dum1(1,1);
m=dum2(1,2);
m_H=dum3(1,2);
q=dum4(1,1);

cons=.001;

d=cons;
g=cons*100;
h=cons;

//beam
M_L=10*cons*[1,0,0,2,1,0]';
N_L=10*cons;
M_K=cons*[1;2];
N_K=cons*[1,0,1,0,1,0];
M_C=cons;
N_H=cons;
M_B=cons*[0,0,3,0,1,0]';
N_B=cons*[1,1];
N_C=cons*[0,2,0,0,1,0];
N1=N_C;
M1=M_B;
N_A=M1';
M_A=M1;
M_H=M1;


save('sys.dat',A,B,C,H,n,m,q,m_H,d,g,h,M_L,N_L,M_K,N_K,M_C,N_C,M_B,N_B,M_H,N_H,M_A,N_A);

function [LME_Hinf, LMI_Hinf, OBJ_Hinf]=nonfragile_Hinf(XLIST_Hinf)
[P,Ph,R,Kh,Lh,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,gama1s]= XLIST_Hinf(:)
LME_Hinf = list(P*B-B*Ph,P-P',R-R',Ph-Ph')
LMI_Hinf = list(-([A'*P+P*A+Kh'*B'+B*Kh+C'*C+(eps4+eps5+eps8+d^2)*eye(n,n),-B*Kh,P*H,2*P*M_B,eps12*N_K',eps1*N_K',P*B*M_K,g*P,h*P,C',P*M_B,zeros(n,1+1+n+n+q+1),P,zeros(n,n);-Kh'*B',A'*R+R*A+C'*Lh'+Lh*C+eps9*d^2*eye(n,n)+eps10*N_C'*N_C,R*H,-2*P*M_B,-eps12*N_K',-eps1*N_K',zeros(n,1+n+n+q+1),R*M_L,eps1*C'*N_L',g*R,h*R,Lh,R*M_L,zeros(n,n),R;H'*P,H'*R,(eps6+eps7-gama1s)*eye(m_H,m_H),zeros(m_H,1+1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);2*M_B'*P,-2*M_B'*P,zeros(1,m_H),-eps11*eye(1,1),zeros(1,1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);eps12*N_K,-eps12*N_K,zeros(1,m_H+1),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);eps1*N_K,-eps1*N_K,zeros(1,m_H+1+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);M_K'*B'*P,zeros(1,n+m_H+1+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);g*P,zeros(n,n+m_H+1+1+1+1),-eps4*eye(n,n),zeros(n,n+q+1+1+1+n+n+q+1),zeros(n,n+n);h*P,zeros(n,n+m_H+1+1+1+1+n),-eps6*eye(n,n),zeros(n,q+1+1+1+n+n+q+1),zeros(n,n+n);C,zeros(q,n+m_H+1+1+1+1+n+n),-eps8*eye(q,q),zeros(q,1+1+1+n+n+q+1),zeros(q,n+n);M_B'*P,zeros(1,n+m_H+1+1+1+1+n+n+q),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1),zeros(1,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1),zeros(1,n+n);zeros(1,n),eps1*N_L*C,zeros(1,m_H+1+1+1+1+n+n+q+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1),zeros(1,n+n);zeros(n,n),g*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1),-eps5*eye(n,n),zeros(n,n+q+1),zeros(n,n+n);zeros(n,n),h*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n),-eps7*eye(n,n),zeros(n,q+1),zeros(n,n+n);zeros(q,n),Lh',zeros(q,m_H+1+1+1+1+n+n+q+1+1+1+n+n),-eps9*eye(q,q),zeros(q,1),zeros(q,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q),-eps10*eye(1,1),zeros(1,n+n);P,zeros(n,n+m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1),-eps2*eye(n,n),zeros(n,n);zeros(n,n),R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1+n),-eps3*eye(n,n)]),P,R,Ph,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,gama1s-0.99)
OBJ_Hinf = gama1s
endfunction

P0=eye(n,n);
Ph0=eye(m,m);
R0=eye(n,n);
Kh0=rand(m,n);
Lh0=rand(n,q);
eps1_0=1;
eps2_0=1;
eps3_0=1;
eps4_0=1;
eps5_0=1;
eps6_0=1;
eps7_0=1;
eps8_0=1;
eps9_0=1;
eps10_0=1;
eps11_0=1;
eps12_0=1;
gama1s0=1;
Init_guess=list(P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,gama1s0);

Ans_LMI_Hinf=lmisolver(Init_guess,nonfragile_Hinf);


P0=Ans_LMI_Hinf(1);
Ph0=Ans_LMI_Hinf(2);
R0=Ans_LMI_Hinf(3);
Kh0=Ans_LMI_Hinf(4);
Lh0=Ans_LMI_Hinf(5);
eps1_0=Ans_LMI_Hinf(6);
eps2_0=Ans_LMI_Hinf(7);
eps3_0=Ans_LMI_Hinf(8);
eps4_0=Ans_LMI_Hinf(9);
eps5_0=Ans_LMI_Hinf(10);
eps6_0=Ans_LMI_Hinf(11);
eps7_0=Ans_LMI_Hinf(12);
eps8_0=Ans_LMI_Hinf(13);
eps9_0=Ans_LMI_Hinf(14);
eps10_0=Ans_LMI_Hinf(15);
eps11_0=Ans_LMI_Hinf(16);
eps12_0=Ans_LMI_Hinf(17);
gama1s0=Ans_LMI_Hinf(18);

K0=Ph0^-1*Kh0;

save('init_vals.dat',P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,gama1s0,K0);

clc
clear
load('init_vals.dat','P0','Ph0','R0','Kh0','Lh0','eps1_0','eps2_0','eps3_0','eps4_0','eps5_0','eps6_0','eps7_0','eps8_0','eps9_0','eps10_0','eps11_0','eps12_0','gama1s0','K0');
load('sys.dat','A','B','C','H','n','m','q','m_H','d','g','h','M_L','N_L','M_K','N_K','M_C','N_C','M_B','N_B','M_H','N_H','M_A','N_A');

C1=zeros(m,n);

dum5=size(C1);
q_C1=dum5(1,1);

function [LME, LMI, OBJ]=nonfragile(XLIST)
[P,Ph,R,Kh,Lh,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,ksi1,ksi2,ksi3,ksi4,ksi5,ksi6,ksi7,ksi8,ksi9,ksi10,ksi11,Z,gama1s,gama2s]= XLIST(:)
LME = list(P*B-B*Ph,P-P',R-R',Ph-Ph',Z-Z')
LMI = list(-([A'*P+P*A+Kh'*B'+B*Kh+C'*C+(eps4+eps5+eps8+d^2)*eye(n,n),-B*Kh,P*H,2*P*M_B,eps12*N_K',eps1*N_K',P*B*M_K,g*P,h*P,C',P*M_B,zeros(n,1+1+n+n+q+1),P,zeros(n,n);-Kh'*B',A'*R+R*A+C'*Lh'+Lh*C+eps9*d^2*eye(n,n)+eps10*N_C'*N_C,R*H,-2*P*M_B,-eps12*N_K',-eps1*N_K',zeros(n,1+n+n+q+1),R*M_L,eps1*C'*N_L',g*R,h*R,Lh,R*M_L,zeros(n,n),R;H'*P,H'*R,(eps6+eps7-gama1s)*eye(m_H,m_H),zeros(m_H,1+1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);2*M_B'*P,-2*M_B'*P,zeros(1,m_H),-eps11*eye(1,1),zeros(1,1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);eps12*N_K,-eps12*N_K,zeros(1,m_H+1),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);eps1*N_K,-eps1*N_K,zeros(1,m_H+1+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);M_K'*B'*P,zeros(1,n+m_H+1+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);g*P,zeros(n,n+m_H+1+1+1+1),-eps4*eye(n,n),zeros(n,n+q+1+1+1+n+n+q+1),zeros(n,n+n);h*P,zeros(n,n+m_H+1+1+1+1+n),-eps6*eye(n,n),zeros(n,q+1+1+1+n+n+q+1),zeros(n,n+n);C,zeros(q,n+m_H+1+1+1+1+n+n),-eps8*eye(q,q),zeros(q,1+1+1+n+n+q+1),zeros(q,n+n);M_B'*P,zeros(1,n+m_H+1+1+1+1+n+n+q),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1),zeros(1,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1),zeros(1,n+n);zeros(1,n),eps1*N_L*C,zeros(1,m_H+1+1+1+1+n+n+q+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1),zeros(1,n+n);zeros(n,n),g*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1),-eps5*eye(n,n),zeros(n,n+q+1),zeros(n,n+n);zeros(n,n),h*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n),-eps7*eye(n,n),zeros(n,q+1),zeros(n,n+n);zeros(q,n),Lh',zeros(q,m_H+1+1+1+1+n+n+q+1+1+1+n+n),-eps9*eye(q,q),zeros(q,1),zeros(q,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q),-eps10*eye(1,1),zeros(1,n+n);P,zeros(n,n+m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1),-eps2*eye(n,n),zeros(n,n);zeros(n,n),R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1+n),-eps3*eye(n,n)]),-([A'*P+P*A+B*Kh+Kh'*B'+ksi1*N_A'*N_A,-B*Kh,C1',ksi2*N_K',ksi4*N_K',ksi9*N_K',P*M_A,P*B*M_K,P*M_B,P*M_B,zeros(n,1+1+1+1+1);-Kh'*B',A'*R+R*A+C'*Lh'+Lh*C+ksi7*C'*N_L'*N_L*C+ksi6*N_C'*N_C+ksi5*N_A'*N_A+ksi8*N_C'*N_C,zeros(n,q_C1),-ksi2*N_K',-ksi4*N_K',-ksi9*N_K',zeros(n,1+1+1+1),R*M_A,Lh*M_C,R*M_L,R*M_L,zeros(n,1);C1,zeros(q_C1,n),-eye(q_C1,q_C1),zeros(q_C1,12);ksi2*N_K,-ksi2*N_K,zeros(1,q_C1),-ksi2*eye(1,1),zeros(1,11);ksi4*N_K,-ksi4*N_K,zeros(1,q_C1+1),-ksi4*eye(1,1),zeros(1,10);ksi9*N_K,-ksi9*N_K,zeros(1,q_C1+2),-ksi9*eye(1,1),zeros(1,9);M_A'*P,zeros(1,n+q_C1+3),-ksi1*eye(1,1),zeros(1,8);M_K'*B'*P,zeros(1,n+q_C1+4),-ksi2*eye(1,1),zeros(1,7);M_B'*P,zeros(1,n+q_C1+5),-ksi3*eye(1,1),zeros(1,6);M_B'*P,zeros(1,n+q_C1+6),-ksi4*eye(1,1),zeros(1,5);zeros(1,n),M_A'*R,zeros(1,q_C1+7),-ksi5*eye(1,1),zeros(1,4);zeros(1,n),M_C'*Lh',zeros(1,q_C1+8),-ksi6*eye(1,1),zeros(1,3);zeros(1,n),M_L'*R,zeros(1,q_C1+9),-ksi7*eye(1,1),zeros(1,2);zeros(1,n),M_L'*R,zeros(1,q_C1+10),-ksi8*eye(1,1),zeros(1,1);zeros(1,n+n+q_C1+11),-ksi9*eye(1,1)]),[P,zeros(n,n),P*H,P*M_H,zeros(n,1);zeros(n,n),R,R*H,zeros(n,1),R*M_H;H'*P,H'*R,Z+ksi11*N_H'*N_H+ksi10*N_H'*N_H,zeros(m_H,2);M_H'*P,zeros(1,n+m_H),-ksi10*eye(1,1),zeros(1,1);zeros(1,n),M_H'*R,zeros(1,m_H+1),-ksi11*eye(1,1)],-(trace(Z)-gama2s),P,R,Ph,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,ksi1,ksi2,ksi3,ksi4,ksi5,ksi6,ksi7,ksi8,ksi9,-ksi10,-ksi11,Z,gama1s,gama2s,gama1s-0.90)
OBJ = gama1s*1
endfunction

Z0=eye(q,q);
ksi1_0=100;
ksi2_0=100;
ksi3_0=100;
ksi4_0=100;
ksi5_0=100;
ksi6_0=100;
ksi7_0=100;
ksi8_0=100;
ksi9_0=100;
ksi10_0=100;
ksi11_0=100;
gama2s0=1;
Init_guess=list(P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,ksi1_0,ksi2_0,ksi3_0,ksi4_0,ksi5_0,ksi6_0,ksi7_0,ksi8_0,ksi9_0,ksi10_0,ksi11_0,Z0,gama1s0,gama2s0);

Ans_LMI=lmisolver(Init_guess,nonfragile);


P0=Ans_LMI(1);
Ph0=Ans_LMI(2);
R0=Ans_LMI(3);
Kh0=Ans_LMI(4);
Lh0=Ans_LMI(5);
eps1_0=Ans_LMI(6);
eps2_0=Ans_LMI(7);
eps3_0=Ans_LMI(8);
eps4_0=Ans_LMI(9);
eps5_0=Ans_LMI(10);
eps6_0=Ans_LMI(11);
eps7_0=Ans_LMI(12);
eps8_0=Ans_LMI(13);
eps9_0=Ans_LMI(14);
eps10_0=Ans_LMI(15);
eps11_0=Ans_LMI(16);
eps12_0=Ans_LMI(17);
ksi1_0=Ans_LMI(18);
ksi2_0=Ans_LMI(19);
ksi3_0=Ans_LMI(20);
ksi4_0=Ans_LMI(21);
ksi5_0=Ans_LMI(22);
ksi6_0=Ans_LMI(23);
ksi7_0=Ans_LMI(24);
ksi8_0=Ans_LMI(25);
ksi9_0=Ans_LMI(26);
ksi10_0=Ans_LMI(27);
ksi11_0=Ans_LMI(28);
Z0=Ans_LMI(29);
gama1s0=Ans_LMI(30);
gama2s0=Ans_LMI(31);

K0=Ph0^-1*Kh0;
L0=R0^-1*Lh0;

save('init_vals1.dat',P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,K0,ksi1_0,ksi2_0,ksi3_0,ksi4_0,ksi5_0,ksi6_0,ksi7_0,ksi8_0,ksi9_0,ksi10_0,ksi11_0,Z0,gama1s0,gama2s0);

endung=10

for ii=1:endung
clc
clear

load('init_vals1.dat','P0','Ph0','R0','Kh0','Lh0','eps1_0','eps2_0','eps3_0','eps4_0','eps5_0','eps6_0','eps7_0','eps8_0','eps9_0','eps10_0','eps11_0','eps12_0','K0','ksi1_0','ksi2_0','ksi3_0','ksi4_0','ksi5_0','ksi6_0','ksi7_0','ksi8_0','ksi9_0','ksi10_0','ksi11_0','Z0','gama1s0','gama2s0');
load('sys.dat','A','B','C','H','n','m','q','m_H','d','g','h','M_L','N_L','M_K','N_K','M_C','N_C','M_B','N_B','M_H','N_H','M_A','N_A');

gama1s=gama1s0;

C1=zeros(m,n);
D1=zeros(m,m);
dum5=size(C1);
q_C1=dum5(1,1);

function [LME, LMI, OBJ]=nonfragile(XLIST)
[P,Ph,R,Kh,Lh,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,ksi1,ksi2,ksi3,ksi4,ksi5,ksi6,ksi7,ksi8,ksi9,ksi10,ksi11,Z,gama2s]= XLIST(:)
LME = list(P*B-B*Ph,P-P',R-R',Ph-Ph',Z-Z')
LMI = list(-([A'*P+P*A+Kh'*B'+B*Kh+C'*C+(eps4+eps5+eps8+d^2)*eye(n,n),-B*Kh,P*H,2*P*M_B,eps12*N_K',eps1*N_K',P*B*M_K,g*P,h*P,C',P*M_B,zeros(n,1+1+n+n+q+1),P,zeros(n,n);-Kh'*B',A'*R+R*A+C'*Lh'+Lh*C+eps9*d^2*eye(n,n)+eps10*N_C'*N_C,R*H,-2*P*M_B,-eps12*N_K',-eps1*N_K',zeros(n,1+n+n+q+1),R*M_L,eps1*C'*N_L',g*R,h*R,Lh,R*M_L,zeros(n,n),R;H'*P,H'*R,(eps6+eps7-gama1s)*eye(m_H,m_H),zeros(m_H,1+1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);2*M_B'*P,-2*M_B'*P,zeros(1,m_H),-eps11*eye(1,1),zeros(1,1+1+1+n+n+q+1+1+1+n+n+q+1),zeros(m_H,n+n);eps12*N_K,-eps12*N_K,zeros(1,m_H+1),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);eps1*N_K,-eps1*N_K,zeros(1,m_H+1+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);M_K'*B'*P,zeros(1,n+m_H+1+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1+1+1+n+n+q+1),zeros(1,n+n);g*P,zeros(n,n+m_H+1+1+1+1),-eps4*eye(n,n),zeros(n,n+q+1+1+1+n+n+q+1),zeros(n,n+n);h*P,zeros(n,n+m_H+1+1+1+1+n),-eps6*eye(n,n),zeros(n,q+1+1+1+n+n+q+1),zeros(n,n+n);C,zeros(q,n+m_H+1+1+1+1+n+n),-eps8*eye(q,q),zeros(q,1+1+1+n+n+q+1),zeros(q,n+n);M_B'*P,zeros(1,n+m_H+1+1+1+1+n+n+q),-eps12*eye(1,1),zeros(1,1+1+n+n+q+1),zeros(1,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1),-eps1*eye(1,1),zeros(1,1+n+n+q+1),zeros(1,n+n);zeros(1,n),eps1*N_L*C,zeros(1,m_H+1+1+1+1+n+n+q+1+1),-eps1*eye(1,1),zeros(1,n+n+q+1),zeros(1,n+n);zeros(n,n),g*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1),-eps5*eye(n,n),zeros(n,n+q+1),zeros(n,n+n);zeros(n,n),h*R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n),-eps7*eye(n,n),zeros(n,q+1),zeros(n,n+n);zeros(q,n),Lh',zeros(q,m_H+1+1+1+1+n+n+q+1+1+1+n+n),-eps9*eye(q,q),zeros(q,1),zeros(q,n+n);zeros(1,n),M_L'*R,zeros(1,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q),-eps10*eye(1,1),zeros(1,n+n);P,zeros(n,n+m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1),-eps2*eye(n,n),zeros(n,n);zeros(n,n),R,zeros(n,m_H+1+1+1+1+n+n+q+1+1+1+n+n+q+1+n),-eps3*eye(n,n)]),-([A'*P+P*A+B*Kh+Kh'*B'+ksi1*N_A'*N_A,-B*Kh,C1'+K0'*D1',ksi2*N_K',ksi4*N_K',ksi9*N_K',P*M_A,P*B*M_K,P*M_B,P*M_B,zeros(n,1+1+1+1+1);-Kh'*B',A'*R+R*A+C'*Lh'+Lh*C+ksi7*C'*N_L'*N_L*C+ksi6*N_C'*N_C+ksi5*N_A'*N_A+ksi8*N_C'*N_C,-K0'*D1',-ksi2*N_K',-ksi4*N_K',-ksi9*N_K',zeros(n,1+1+1+1),R*M_A,Lh*M_C,R*M_L,R*M_L,zeros(n,1);C1+D1*K0,-D1*K0,-eye(q_C1,q_C1),zeros(q_C1,11),D1*M_K;ksi2*N_K,-ksi2*N_K,zeros(1,q_C1),-ksi2*eye(1,1),zeros(1,11);ksi4*N_K,-ksi4*N_K,zeros(1,q_C1+1),-ksi4*eye(1,1),zeros(1,10);ksi9*N_K,-ksi9*N_K,zeros(1,q_C1+2),-ksi9*eye(1,1),zeros(1,9);M_A'*P,zeros(1,n+q_C1+3),-ksi1*eye(1,1),zeros(1,8);M_K'*B'*P,zeros(1,n+q_C1+4),-ksi2*eye(1,1),zeros(1,7);M_B'*P,zeros(1,n+q_C1+5),-ksi3*eye(1,1),zeros(1,6);M_B'*P,zeros(1,n+q_C1+6),-ksi4*eye(1,1),zeros(1,5);zeros(1,n),M_A'*R,zeros(1,q_C1+7),-ksi5*eye(1,1),zeros(1,4);zeros(1,n),M_C'*Lh',zeros(1,q_C1+8),-ksi6*eye(1,1),zeros(1,3);zeros(1,n),M_L'*R,zeros(1,q_C1+9),-ksi7*eye(1,1),zeros(1,2);zeros(1,n),M_L'*R,zeros(1,q_C1+10),-ksi8*eye(1,1),zeros(1,1);zeros(1,n+n),M_K'*D1',zeros(1,11),-ksi9*eye(1,1)]),[P,zeros(n,n),P*H,P*M_H,zeros(n,1);zeros(n,n),R,R*H,zeros(n,1),R*M_H;H'*P,H'*R,Z+ksi11*N_H'*N_H+ksi10*N_H'*N_H,zeros(m_H,2);M_H'*P,zeros(1,n+m_H),-ksi10*eye(1,1),zeros(1,1);zeros(1,n),M_H'*R,zeros(1,m_H+1),-ksi11*eye(1,1)],-(trace(Z)-gama2s),P,R,Ph,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,eps11,eps12,ksi1,ksi2,ksi3,ksi4,ksi5,ksi6,ksi7,ksi8,ksi9,-ksi10,-ksi11,Z,gama2s)
OBJ = []
endfunction

Init_guess=list(P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,ksi1_0,ksi2_0,ksi3_0,ksi4_0,ksi5_0,ksi6_0,ksi7_0,ksi8_0,ksi9_0,ksi10_0,ksi11_0,Z0,gama2s0);

Ans_LMI=lmisolver(Init_guess,nonfragile);


P0=Ans_LMI(1);
Ph0=Ans_LMI(2);
R0=Ans_LMI(3);
Kh0=Ans_LMI(4);
Lh0=Ans_LMI(5);
eps1_0=Ans_LMI(6);
eps2_0=Ans_LMI(7);
eps3_0=Ans_LMI(8);
eps4_0=Ans_LMI(9);
eps5_0=Ans_LMI(10);
eps6_0=Ans_LMI(11);
eps7_0=Ans_LMI(12);
eps8_0=Ans_LMI(13);
eps9_0=Ans_LMI(14);
eps10_0=Ans_LMI(15);
eps11_0=Ans_LMI(16);
eps12_0=Ans_LMI(17);
ksi1_0=Ans_LMI(18);
ksi2_0=Ans_LMI(19);
ksi3_0=Ans_LMI(20);
ksi4_0=Ans_LMI(21);
ksi5_0=Ans_LMI(22);
ksi6_0=Ans_LMI(23);
ksi7_0=Ans_LMI(24);
ksi8_0=Ans_LMI(25);
ksi9_0=Ans_LMI(26);
ksi10_0=Ans_LMI(27);
ksi11_0=Ans_LMI(28);
Z0=Ans_LMI(29);
gama2s0=Ans_LMI(30);

K_ind=K0-Ph0^-1*Kh0;
K0=Ph0^-1*Kh0;
L0=R0^-1*Lh0;

save('init_vals1.dat',P0,Ph0,R0,Kh0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,eps11_0,eps12_0,K0,ksi1_0,ksi2_0,ksi3_0,ksi4_0,ksi5_0,ksi6_0,ksi7_0,ksi8_0,ksi9_0,ksi10_0,ksi11_0,Z0,gama1s0,gama2s0);
end
K1=K0;
L1=L0;
save('L1.dat',L1);
save('K1.dat',K1);

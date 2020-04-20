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

//beam
M_L=10*cons*[1,0,0,2,1,0]';
N_L=10*cons;
M_K=cons*[1;2];
N_K=cons*[1,0,1,0,1,0];

row=.48640;


function [LME_Hinf, LMI_Hinf, OBJ_Hinf]=nonfragile_Hinf(XLIST_Hinf)
[P,Ph,R,Kh,Lh,eps1]= XLIST_Hinf(:)
LME_Hinf = list(P*B-B*Ph,P-P',R-R',Ph-Ph')
LMI_Hinf = list(-([A'*P+P*A-Kh'*B'-B*Kh+2*row*P,B*Kh,P*B*M_K,zeros(n,1),-eps1*N_K',zeros(n,1);Kh'*B',A'*R+R*A-Lh*C-C'*Lh'+2*row*R,zeros(n,1),R*M_L,eps1*N_K',-eps1*C'*N_L';M_K'*B'*P',zeros(1,n),-eps1,zeros(1,1+1+1);zeros(1,n),M_L'*R,0,-eps1,zeros(1,2);-eps1*N_K,eps1*N_K,zeros(1,2),-eps1,0;zeros(1,n),-eps1*N_L*C,zeros(1,3),-eps1]),P,R,Ph,eps1)
OBJ_Hinf = []
endfunction

P0=eye(n,n);
Ph0=eye(m,m);
R0=eye(n,n);
Kh0=rand(m,n);
Lh0=rand(n,q);
eps1_0=1;

Init_guess=list(P0,Ph0,R0,Kh0,Lh0,eps1_0);

Ans_LMI_Hinf=lmisolver(Init_guess,nonfragile_Hinf);


P=Ans_LMI_Hinf(1);
Ph=Ans_LMI_Hinf(2);
R=Ans_LMI_Hinf(3);
Kh=Ans_LMI_Hinf(4);
Lh=Ans_LMI_Hinf(5);
eps1=Ans_LMI_Hinf(6);

K=Ph^-1*Kh;
L=R^-1*Lh;

K2=K;
L2=L;
save('L2.dat',L2);
save('K2.dat',K2);

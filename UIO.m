clc;
clear all
close all
warning('off')
format shortE
%%
AC=[-0.0011, 0.0769,0,0;
    0       -2.5, 2.5,0;
    -5.2083,0,-12.5,12.5;
    -0.301 0 0 0];

BC=[-0.0769;
    0;
    0;
    0];

C=[1 0 0 0
    0 0 0 1];

D=[0
   0];


sys_ss = ss(AC,BC,C,D);
Ts = 2;
sys_d=c2d(sys_ss,Ts,'zoh');
A=sys_d.a
B=sys_d.b

x=size(A,1);
m=size(B,2);
p=size(C,1);

M1=[D,zeros(size(D));C*B,D];
M2=M1;
flag=1;
theta_alpha_1=[C];
R=-1;

while R~=m
    M1=M2;
    theta_alpha_1=[theta_alpha_1;C*A^(flag)];
    M2=[D,zeros(size(D,1),size(M1,2));theta_alpha_1*B,M1];
    flag=flag+1;
    R=rank(M2)-rank(M1);
end

L=flag;
ML=M2;

N1=left_null_space(M1');
N2=[eye(p),zeros(p,size(N1,2));zeros(size(N1,1),p),N1];
W=[zeros(m,m);eye(m)]*pinv(N2*[D;theta_alpha_1*B]);

WW=left_null_space([N2*[D;theta_alpha_1*B]]');
W(1:m,:)=WW(1:m,:);

Q=W*N2
Q*M2;

theta_alpha=[theta_alpha_1;C*A^L];
S=Q*theta_alpha;
SS=size(S,1);
S1=S(1:SS/2,:)
S2=S(SS/2+1:end,:)

New_Pole_Places=[-.1,-0.05,0,0.01];
L1=place((A-B*S2)',S1',New_Pole_Places)
L2=B
L=[L1',L2]*Q




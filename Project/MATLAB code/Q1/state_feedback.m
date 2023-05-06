clc
clear
studentID = 'A0000001X';
parameters = get_parameter();
A = parameters{1};
B = parameters{2};
C = parameters{3};
x0 = parameters{5};

%% reference model
syms s;
ts=5;%settling time
mp=0.1;%overshoot

ep_min=abs(log(mp)/sqrt(pi^2+(log(mp))^2));
ep=0.8;
wn=4/(ep*ts);

% lamda1=-ep*wn+wn*sqrt(1-ep^2)*1i;
% lamda2=-ep*wn-wn*sqrt(1-ep^2)*1i;
lamda1=-0.8+0.3i;
lamda2=-0.8-0.3i;
lamda3=- 3;
lamda4=- 3;
lamda5=- 2;
lamda6=- 2;

polynomial1=(s-lamda1)*(s-lamda2);
polynomial2=(s-lamda1)*(s-lamda2)*(s-lamda3)*(s-lamda4)*(s-lamda5)*(s-lamda6);
polynomial3 = (s-lamda1)*(s-lamda2)*(s-lamda3);
polynomial4 = (s-lamda4)*(s-lamda5)*(s-lamda6);

Ad_cof1=double(coeffs(polynomial1));
Ad_cof2=double(coeffs(polynomial2));
Ad_cof3=double(coeffs(polynomial3));
Ad_cof4=double(coeffs(polynomial4));
K_verify = place(A,B,[lamda1 lamda2 lamda3 lamda4 lamda5 lamda6]);

%% Full Rank pole placement

syms k11 k12 k13 k14 k15 k16 k21 k22 k23 k24 k25 k26 ;
K=[k11 k12 k13 k14 k15 k16; 
      k21 k22 k23 k24 k25 k26];

% controllability matrix
W=[B A*B A^2*B A^3*B];
% verify if W is controlible
assert(rank(W(:,1:6))==6);

Cmatirx = [B(:,1),A*B(:,1),A^2*B(:,1),B(:,2),A*B(:,2),A^2*B(:,2)];
inv_Cmatrix=inv(Cmatirx);
bcol1=3;
bcol2=3;

T=[inv_Cmatrix(bcol1,:); 
      inv_Cmatrix(bcol1,:)*A;
      inv_Cmatrix(bcol1,:)*A^2;
      inv_Cmatrix(bcol1+bcol2,:);
      inv_Cmatrix(bcol1+bcol2,:)*A;
      inv_Cmatrix(bcol1+bcol2,:)*A^2];

A_bar=T*A/(T);
B_bar=T*B;
A_bar(abs(A_bar)<10^(-10))=0;
B_bar(abs(B_bar)<10^(-10))=0;

Am=A_bar-B_bar*K;
Ad=[0 1 0 0 0 0 ; 
        0 0 1 0 0 0;
        -Ad_cof3(1:3) 0 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        0 0 0 -Ad_cof4(1:3)];

 %% solve rotation
rotation=Am==Ad;
K_num=solve(rotation);
K_ans=struct2array(K_num);
K_ans=double(K_ans);
K_ba=[K_ans(1:6);
    K_ans(7:12)];
K_calculated=K_ba*T;

%% Plot the performence

t=0:0.02:10;
Af=A-B*K_calculated;
sys=ss(Af,B,C,0);

len = size(t,2);
u0=zeros(len,2);
u1=[ones(len,1),zeros(len,1)];
u2=[zeros(len,1),ones(len,1)];

% zero inputs and x0 initial state
[y,tout,x]=lsim(sys,u0,t,x0);

figure()
plot(t,x)
grid on
legend('x1','x2','x3','x4','x5','x6')
xlabel('time')
ylabel('state')
title('zero inputs and x0 initial state')

figure()
plot(t,y)
grid on
legend('y1','y2','y3')
xlabel('time')
ylabel('output')
title('zero inputs and x0 initial state')

for i = 1:length(t)
% calculate the u
u_0(i,:) = -K_calculated*x(i,:)';
end
figure()
plot(t,u_0)
grid on
legend('uc','uh')
xlabel('time')
ylabel('control signal')
title('zero inputs and x0 initial state')

figure()
step(sys);
grid on


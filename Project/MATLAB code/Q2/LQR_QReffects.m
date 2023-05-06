clc;
clear;
% get parameters
parameters = get_parameter();
A = parameters{1};
B = parameters{2};
C = parameters{3};
C2 = parameters{4};
x0 = parameters{5};
D = parameters{6};
ysp = parameters{7};

%% evaluate LQR performance when changing Q and R 
Q1=[1 0 0 0 0 0
   0 1 0 0 0 0
   0 0 1 0 0 0
   0 0 0 1 0 0
   0 0 0 0 1 0
   0 0 0 0 0 1].*50;
%R=[1 0
   %0 1].*0.5;

Q2=[1 0 0 0 0 0
   0 1 0 0 0 0
   0 0 1 0 0 0
   0 0 0 1 0 0
   0 0 0 0 1 0
   0 0 0 0 0 1].*0.5;
R=[1 0
   0 1];

% change Q
gamma1=[A -B/R*B';-Q1 -A'];
gamma2=[A -B/R*B';-Q2 -A'];

[eig_vector1,eig_value1]=eig(gamma1);
eig_value1_sum=sum(eig_value1);
vueigen1=eig_vector1(:,real(eig_value1_sum)<0);
P1=vueigen1(7:12,:)/vueigen1(1:6,:);
K1=inv(R)*B'*P1;

[eig_vector2,eig_value2]=eig(gamma2);
eig_value2_sum=sum(eig_value2);
vueigen2=eig_vector2(:,real(eig_value2_sum)<0);
P2=vueigen2(7:12,:)/vueigen2(1:6,:);
K2=inv(R)*B'*P2;

%% PLOT figure
t=0:0.01:10;
Af1=A-B*K1;
Af2=A-B*K2;
sys1=ss(Af1,B,C,D);
sys2=ss(Af2,B,C,D);

len = size(t,2);
u0=zeros(len,2);
u1=[ones(len,1),zeros(len,1)];
u2=[zeros(len,1),ones(len,1)];

[y1,tout1,x1]=lsim(sys1,u0,t,x0);
figure()
plot(t,x1)
legend('x1 Q=50I','x2 Q=50I','x3 Q=50I','x4 Q=50I')
xlabel('time')
ylabel('state')
grid on
figure()
plot(t,y1)
legend('out1 Q=50I','out2')
grid on

%hold on
[y2,tout2,x2]=lsim(sys2,u0,t,x0);
figure()
plot(t,x2)
legend('x1 Q=0.5I','x2 Q=0.5I','x3 Q=0.5I','x4 Q=0.5I')
xlabel('time')
ylabel('state')
grid on
figure()
plot(t,y2)
legend('out1 Q=0.5I','out2')
grid on

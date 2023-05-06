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

%% system design

% verify controbility
W=[B A*B A^2*B A^3*B];
assert(rank(W(:,1:6))==6);

% LQR method
Q=[1 0 0 0 0 0
   0 10 0 0 0 0
   0 0 100 0 0 0
   0 0 0 100 0 0
   0 0 0 0 10 0
   0 0 0 0 0 1]*1;
R=[100 0
   0 1];

%[K1,~,P]=lqr(A,B,Q,R)

gamma=[A -B/R*B';-Q -A'];
[eig_vector,eig_value]=eig(gamma);
eig_value_sum=sum(eig_value);
vueigen=eig_vector(:,real(eig_value_sum)<0);
P=vueigen(7:12,:)/vueigen(1:6,:);
K_calculated=real(inv(R)*B'*P);

%% PLOT figure
t=0:0.01:10;
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
legend('x1','x2','x3','x4','x5','x6')
xlabel('time')
ylabel('state')
title('zero inputs and x0 initial state')

figure()
plot(t,y)
legend('y1','y2','y3')
xlabel('time')
ylabel('output')
title('zero inputs and x0 initial state')

for i = 1:length(t)
    u_in(i,:) = -K_calculated*x(i,:)';
end
figure()
plot(t,u_in)
legend('uc','uh')
xlabel('time')
ylabel('control signal')
title('zero inputs and x0 initial state')


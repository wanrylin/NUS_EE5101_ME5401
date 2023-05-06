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

%% Servo control + LQR
w = [-1;1];     

% verify controlbility
Qc = [A B;C zeros(3,2)];
assert(rank(Qc)==8);

A_bar=[A zeros(6,3);-C zeros(3,3)];
B_bar=[B;zeros(3,2)];
B_w_bar=[B;zeros(3,2)];
B_r_bar=[zeros(6,3);eye(3)];
C_bar=[C,zeros(3,3)];

Q=[1 0 0 0 0 0 0 0 0;
   0 1 0 0 0 0 0 0 0;
   0 0 1 0 0 0 0 0 0;
   0 0 0 1 0 0 0 0 0;
   0 0 0 0 1 0 0 0 0;
   0 0 0 0 0 1 0 0 0;
   0 0 0 0 0 0 1 0 0;
   0 0 0 0 0 0 0 1 0;
   0 0 0 0 0 0 0 0 1]*10;
R=[1 0
   0 1]*1;

gamma=[A_bar -B_bar/R*(B_bar');-Q -A_bar'];

[eig_vector,eig_value]=eig(gamma);
eig_value_sum=sum(eig_value);
vueogen=eig_vector(:,real(eig_value_sum)<0);
P=vueogen(10:18,:)/vueogen(1:9,:);
K_calculated=real(inv(R)*(B_bar')*P);

K1=K_calculated(:,1:6);
K2=K_calculated(:,7:9);

a11 = A-B*K1;
a12 = -B*K2;
Au = [a11,a12;-C,zeros(3,3)];
eigvalue = eig(Au);

for i=1:9
np11 = eigvalue(i).*eye(6,6)-A;
np = [np11,B;-C,zeros(3,2)];
a(i)=rank(np);
end

%% full order Observer LQR method
% A_ba = A', B_ba = C', K_ba = L'
Qbar=[1 0 0 0 0 0
      0 1 0 0 0 0
      0 0 1 0 0 0
      0 0 0 1 0 0
      0 0 0 0 1 0
      0 0 0 0 0 1]*5;
Rbar=[1 0 0
       0 1 0
       0 0 1]*1;

Phi1 = [A',-C'/Rbar*C;-Qbar,-A];
[eig_vector_observed,eig_value_observed]=eig(Phi1);
eig_value_observed_sum=sum(eig_value_observed);
vueigen_observered=eig_vector_observed(:,real(eig_value_observed_sum)<0);
P_observed=vueigen_observered(7:12,:)/vueigen_observered(1:6,:);
Kbar=real(inv(Rbar)*C*P_observed);
L=Kbar';
L=real(L);
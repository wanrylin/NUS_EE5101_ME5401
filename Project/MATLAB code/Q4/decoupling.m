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

%% decoupling performance

syms s

for i = 1:6
    if C2(1,:)*A^(i-1)*B ~= 0
        degree1 = i;
        break
    end
end
for i = 1:6
    if C2(2,:)*A^(i-1)*B ~= 0
        degree2 = i;
        break
    end
end

B_star = [C2(1,:)*A^(degree1-1)*B;C2(2,:)*A^(degree2-1)*B];

Phi_A1=A^2+22*A+40*eye(6);

Phi_A2=A^2+50*A+49*eye(6);

C_star = [C2(1,:)*Phi_A1;C2(2,:)*Phi_A2];

F = inv(B_star);
K = F*C_star;
Bf=B*F;
Af = A-B*K;
decouple_model=ss(Af,Bf,C2,0);
W=[Bf Af*Bf Af^2*Bf Af^3*Bf];
assert(rank(W)==6);

p=pole(decouple_model);
H=C2*inv(s*eye(6)-Af)*Bf;

[~,eigenvalue] = eig(Af);

%% Plot
% non-zero inputs and zero initial state
figure()
step(decouple_model);
grid on

% zero inputs and x0 initial state
t=0:0.01:10;
len = size(t,2);
u0=zeros(len,2);
u1=[ones(len,1),zeros(len,1)];
u2=[zeros(len,1),ones(len,1)];

[y,tout,x]=lsim(decouple_model,u0,t,x0);

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
legend('y1','y2')
xlabel('time')
ylabel('output')
title('zero inputs and x0 initial state')


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

%% resultant observer-based LQR system

Q=[1 0 0 0 0 0
   0 1 0 0 0 0
   0 0 1 0 0 0
   0 0 0 1 0 0
   0 0 0 0 1 0
   0 0 0 0 0 1]*10;
R=[1 0
   0 1];

gamma=[A -B/R*B';-Q -A'];
[eig_vector,eig_value]=eig(gamma);
eig_value_sum=sum(eig_value);
vueigen=eig_vector(:,real(eig_value_sum)<0);
P=vueigen(7:12,:)/vueigen(1:6,:);
K_calculated=real(inv(R)*B'*P);
Af=A-B*K_calculated;
[eig_vector_Af,eig_value_Af]=eig(Af);

sys=ss(Af,B,C,D);
orig_pole=pole(sys);

%% full order pole placement,
% A_ba = A', B_ba = C', K_ba = L'
syms s 
dir_pole = [-8,-8,-12,-12,-32,-32];       
polynomial1 = (s-dir_pole(1))*(s-dir_pole(2));
polynomial2 = (s-dir_pole(3))*(s-dir_pole(4));
polynomial3 = (s-dir_pole(5))*(s-dir_pole(6));
Ad_cof1=double(coeffs(polynomial1));
Ad_cof2=double(coeffs(polynomial2));
Ad_cof3=double(coeffs(polynomial3));

syms l11 l12 l13 l14 l15 l16 l21 l22 l23 l24 l25 l26 l31 l32 l33 l34 l35 l36 ;
L=[l11 l12 l13 l14 l15 l16;
    l21 l22 l23 l24 l25 l26;
    l31 l32 l33 l34 l35 l36];

W_bar=[C' A'*C' (A')^2*C' (A')^3*C'];
assert(rank(W_bar(:,1:6))==6);
Cmatrix=W_bar(:,1:6);
C_reconstructed=[Cmatrix(:,1),Cmatrix(:,4),Cmatrix(:,2),Cmatrix(:,5),Cmatrix(:,3),Cmatrix(:,6)];
inv_Cmatrix=inv(C_reconstructed);
bcol1=2;
bcol2=2;
bcol3=2;

T=[inv_Cmatrix(bcol1,:); 
    inv_Cmatrix(bcol1,:)*A';
    inv_Cmatrix(bcol1+bcol2,:);
    inv_Cmatrix(bcol1+bcol2,:)*A';
    inv_Cmatrix(bcol1+bcol2+bcol3,:);
    inv_Cmatrix(bcol1+bcol2+bcol3,:)*A';];

Abar=T*(A')/(T);
Bbar=T*C';
Abar(abs(Abar)<10e-10)=0;
Bbar(abs(Bbar)<10e-10)=0;
Am=Abar-Bbar*L;
Ad=[0 1 0 0 0 0 ; 
        -Ad_cof1(1:2) 0 0 0 0;
        0 0 0 1 0 0;
        0 0 -Ad_cof2(1:2) 0 0;
        0 0 0 0 0 1;
        0 0 0 0 -Ad_cof3(1:2)];

% solve rotation
rotation=Am==Ad;
L_number=solve(rotation);
L_answer=struct2array(L_number);
L_answer=double(L_answer);
Lbar=[L_answer(1:6);
    L_answer(7:12);
    L_answer(13:18)];
L_estimated=Lbar*T;
L_estimated=L_estimated';
L_estimated=real(L_estimated);
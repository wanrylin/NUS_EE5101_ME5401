clc
clear
studentID = 'A0000001X';
parameters = get_parameter();

syms s;
ts= 5 ;%settling time
mp=0.1;%overshoot

ep_min=abs(log(mp)/sqrt(pi^2+(log(mp))^2));
ep=0.8;
wn=4/(ep*ts);

lamda1=-ep*wn+wn*sqrt(1-ep^2)*1i*2;
lamda2=-ep*wn-wn*sqrt(1-ep^2)*1i*2;
lamda3=- 4;
lamda4=- 4;
lamda5=- 6;
lamda6=- 6;

polynomial1=1/(lamda1*lamda2)*(s-lamda1)*(s-lamda2);
polynomial2=1/(lamda1*lamda2*lamda3*lamda4*lamda5*lamda6)*(s-lamda1)*(s-lamda2)*(s-lamda3)*(s-lamda4)*(s-lamda5)*(s-lamda6);
Ad_cof=double(coeffs(polynomial1));
Ad_cof2=double(coeffs(polynomial2));

figure();
step(tf(1,fliplr(Ad_cof2)))
hold on 
step(tf(1,fliplr(Ad_cof)))
legend('add extra pole','original system')
xlim([0,10]);


function output = get_parameter(studentID)
% get studentID: string
% output parameter matrix

    % default studentID is mine
    if nargin < 1
        studentID = 'A0000001X';
    end
    
    % matriculation number
    a = str2num(studentID(5));
    b = str2num(studentID(6));
    c = str2num(studentID(7)); 
    d = str2num(studentID(8));
    
    % Physical parameters
    Mf = 2.14+c/20; 
    Mr = 5.91-b/10;
    Mc = 1.74;
    LFf = 0.05;
    Lr = 0.128; 
    Lc = 0.259;
    Jx = 0.5+(c-d)/100;
    alpha = 15.5-a/3+b/2;
    gamma = 11.5+(a-c)/(b+d+3);

    Hf = 0.18;
    Hr = 0.161;
    Hc = 0.098;
    LF = 0.133;
    LR = 0.308+(a-d)/100;
    mux = 3.33-b/20+a*c/60;
    beta = 27.5-d/2;
    delta = 60+(a-b)*c/10;
    
    den = Mf*Hf*Hf+Mr*Hr*Hr+Mc*Hc*Hc+Jx;
    g = 9.8;
    a51 = -1*Mc*g/den;
    a52 = (Mf*Hf+Mr*Hr+Mc*Hc)*g/den;
    a53 = (Mr*Lr*LF+Mc*Lc*LF+Mf*LFf*LR)*g/((LR+LF)*den);
    a54 = -1*Mc*Hc*alpha/den;
    a55 = -1*mux/den;
    a56 = Mf*Hf*LFf*gamma/den;
    b51 = Mc*Hc*beta/den;
    b52 = -1*Mf*Hf*LFf*delta/den;
    
    A=[0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        0 6.5 -10 -1*alpha 0 0;
        a51 a52 a53 a54 a55 a56;
        5 -3.6 0 0 0 -1*gamma];
    
    B=[0 0;
        0 0;
        0 0;
        beta 11.2;
        b51 b52;
        40 delta];
    
    C=[1 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0];
    
    C2 = [1 0 0 0 0 0;
        0 0 1 0 0 0];
    
    x0=[0.2;-0.1;0.15;-1;0.8;0];
    
    D=[0 0;0 0;0 0];
    
    oprator = [-0.5+(a-b)/20;
        0.1+(b-c)/(a+d+10)];
    y_sp = -0.1*C*inv(A)*B*oprator;

    output = {A,B,C,C2,x0,D,y_sp};

end
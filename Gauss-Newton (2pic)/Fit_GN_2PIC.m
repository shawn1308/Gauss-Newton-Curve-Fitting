clear;clc;
load('Experi.mat');
%% Global VAR
global Ythe; global Xthe; global Yexp; global Xexp; 
global N_DEP; global N_FIN;global Point0;

%% Fitting RANGE - ZEROS
N_DEP = 1;
N_FIN = 1100;

Point0 = 1696;
ec = 1;
%% import
XYthe = importdata('A6B.dat');
Ythe = XYthe(:,2);
Xthe = XYthe(:,1)-Point0;

XYexp = experi;
Yexp = XYexp(:,2);
Xexp = XYexp(:,1);

%% Initial Conditions
Diff_Min_Max_the = max(Ythe) - min(Ythe);
Diff_Min_Max_exp = max(Yexp) - min(Yexp);
[MYthe,Ithe] = min(Ythe); MXthe = Xthe(Ithe);
[MYexp,Iexp] = min(Yexp); MXexp = Xexp(Iexp);

C1 = Diff_Min_Max_exp/Diff_Min_Max_the;% dilat y
C2 = 1; % dilat x
C3 = 10; % shift x
C4 = 0.00001; % offset y
C5 = 0.6e-4;
C6 = 390;

%% AFF sans itération
[fy,fx] = GN_Function_2PIC(C1,C2,C3,C4,C5,C6);
figure;
plot(fx,fy,'g',Xexp(N_DEP:N_FIN),Yexp(N_DEP:N_FIN),'b');
hold on

deltac = 0.1;
deltac2 = 0.0001;

for i=1:100
    fy = GN_Function_2PIC(C1,C2,C3,C4,C5,C6);
    j(:,1) = (GN_Function_2PIC(C1+deltac,C2,C3,C4,C5,C6) - fy)/deltac;
    j(:,2) = (GN_Function_2PIC(C1,C2+deltac2,C3,C4,C5,C6) - fy)/deltac2;
    j(:,3) = (GN_Function_2PIC(C1,C2,C3+deltac,C4,C5,C6) - fy)/deltac;
    j(:,4) = (GN_Function_2PIC(C1,C2,C3,C4+deltac,C5,C6) - fy)/deltac;
    j(:,5) = (GN_Function_2PIC(C1,C2,C3,C4,C5+deltac,C6) - fy)/deltac;
    j(:,6) = (GN_Function_2PIC(C1,C2,C3,C4,C5,C6+deltac) - fy)/deltac;
    d=Yexp(N_DEP:N_FIN)-fy;
    dp = ((j'*j)^(-1))*(j'*d);
    C1 = C1 + dp(1);
    C2 = C2 + dp(2);
    C3 = C3 + dp(3);
    C4 = C4 + dp(4);
    C5 = C5 + dp(5);
    C6 = C6 + dp(6);
    if(abs(dp(1))< 0.01 && abs(dp(2))< 0.01 && abs(dp(3))< 0.01 && abs(dp(4))< 0.01)
        break;
    end
end

%% Affichage Fit
[fy,fx] = GN_Function_2PIC(C1,C2,C3,C4,C5,C6);
plot(fx,fy,'r');

%% Params 
gamma = 2.0*ec*C2/0.1
shift1 = C3
shift2 = C6

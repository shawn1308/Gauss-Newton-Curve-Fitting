clear;clc;
load('Amalan.mat');
%% Global VAR
global Ythe; global Xthe; global Yexp; global Xexp; 
global N_DEP; global N_FIN;global Point0;

%% Fitting RANGE - ZEROS
N_DEP = 900; % intervalle [N_DEP:N_FIN] à fitter (points)
N_FIN = 1200;

Point0 = 1696; % Zero de la courbe theorique
ec = 1;
%% Import DATA
XYthe = importdata('A6B.dat');
Ythe = XYthe(:,2);
Xthe = XYthe(:,1)-Point0; % centrer au point Zero

XYexp = exp;
Yexp = XYexp(:,2);
Xexp = XYexp(:,1);

%% Axe delta
Dlt_the = (max(Xthe) - min(Xthe))/length(Xthe); % (non-utilisé)-Delta fréquence
Dlt_exp = (max(Xexp) - min(Xexp))/length(Xexp);

%% Initial Conditions
Diff_Min_Max_the = max(Ythe) - min(Ythe);
Diff_Min_Max_exp = max(Yexp) - min(Yexp);
[MYthe,Ithe] = min(Ythe); MXthe = Xthe(Ithe);
[MYexp,Iexp] = min(Yexp); MXexp = Xexp(Iexp);

C1 = Diff_Min_Max_exp/Diff_Min_Max_the;% Amplitude coef
C2 = 1; % élargissement x
C3 = 1; % shift x
C4 = mean(Yexp(1:20))/2; % offset y

%% Affichage sans iteration
[fy,fx] = GN_Function_P0(C1,C2,C3,C4);
figure;
plot(fx,fy,'g',Xexp(N_DEP:N_FIN),Yexp(N_DEP:N_FIN),'b');
hold on

%% GN iterative
deltac = 0.01;
for i=1:100
    fy = GN_Function_P0(C1,C2,C3,C4);
    j(:,1) = (GN_Function_P0(C1+deltac,C2,C3,C4) - fy)/deltac; % Jacobien = derivé partielle/parametres
    j(:,2) = (GN_Function_P0(C1,C2+deltac,C3,C4) - fy)/deltac;
    j(:,3) = (GN_Function_P0(C1,C2,C3+deltac,C4) - fy)/deltac;
    j(:,4) = (GN_Function_P0(C1,C2,C3,C4+deltac) - fy)/deltac;
    d=Yexp(N_DEP:N_FIN)-fy;
    dp = ((j'*j)^(-1))*(j'*d);
    C1 = C1 + dp(1);
    C2 = C2 + dp(2);
    C3 = C3 + dp(3);
    C4 = C4 + dp(4);
    if(abs(dp(1))< 0.01 && abs(dp(2))< 0.01 && abs(dp(3))< 0.01 && abs(dp(4))< 0.01) %condition d'arret
        break;
    end
end

%% Affichage Fit
[fy,fx] = GN_Function_P0(C1,C2,C3,C4);
plot(fx,fy,'r');

%% Params
gamma = 2.0*ec*C2/0.1
shift = C3

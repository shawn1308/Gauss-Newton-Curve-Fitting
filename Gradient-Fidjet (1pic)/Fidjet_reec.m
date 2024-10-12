%% Nombre de points
Nthe = 3451; %Nombre de points theorique
Nexp = 2181; %Nombre de points exprimentaux

%% Intervall sur la courbe théorique à Fitter
Ndep = 900; % point de Départ dans l'array expérimentale
Nfin = 1300; % point de Fin

%% Echelle
Dexp = 1; % Intervalle entre 2 point expérimantaux (MHz)
Point0 = 1691; % Zero de courbe théorique
Position0 = 1121; % ~Zéro de courbe expérimentale

%% Import courbes
Data_theorique = importdata('A6B.dat');
Ythe = Data_theorique(:,2);

Data_experimental = exp(:,2).';
Yexp = Data_experimental;

%% Valeur Initiaux - calculé en utilisant min et max
[Ythe_Max,Xthe_Max] = max(Ythe);
[Ythe_Min,Xthe_Min] = min(Ythe);

[Yexp_Max,Xexp_Max] = max(Yexp);
[Yexp_Min,Xexp_Min] = min(Yexp);

C1 = (Yexp_Max - Yexp_Min)/(Ythe_Max - Ythe_Min);
C2 = (Xexp_Max - Xexp_Min)/(Xthe_Max - Xthe_Min);
C3 = (Xexp_Max-Xthe_Max);
C4 = 0;

%% itération
Niter = 5;
dp = zeros(4,1);
for iteration = 1:Niter
    %% Deriver premier et seconde de courbe theorique
    YPthe = diff(Ythe);
    YSthe = diff(YPthe);
    
    %% Linéarisation de fonction par extrapolation
    for Ifunc = Ndep:Nfin
        xtil = (double(Ifunc)-C3)/C2;
        n = int16(xtil);
        d = xtil - double(n);
        if (n+1 <= Nthe) && (n >= 1)
            Yfun(Ifunc) = (Ythe(n+1) - Ythe(n))*d + Ythe(n);
            YPfun(Ifunc) = (YPthe(n+1) - YPthe(n))*d + YPthe(n);
            YSfun(Ifunc) = (YSthe(n+1) - YSthe(n))*d + YSthe(n);
        end
        %% Dérivé partielle par rappor à C2
        DP_C2(Ifunc) = -YPfun(Ifunc)*xtil/C2;
        DS_C2(Ifunc) = YSfun(Ifunc)*xtil.^2/(C2*C2) + 2*YSfun(Ifunc)/(C2*C2)*xtil;
    end

    %% Initialisation C4 par moyennage
    if iteration == 1
        somme_C4 = sum(Yexp(Ndep:Nfin)-C1*Yfun(Ndep:Nfin));
        C4 = somme_C4/(Nfin-Ndep+1);
    end

    %% Calcul erreur quadratique
    

    %% Rempli matrice demi
    A = zeros(4,5);
    for Imatrix = Ndep:Nfin
        AA1 = Yexp(Imatrix)-C1*Yfun(Imatrix)-C4;
        CO1 = double(Imatrix)-C3;
        
        %% Ligne 1
        A(1,1) = A(1,1) + 2.0*Yfun(Imatrix)*Yfun(Imatrix);
        A(1,2) = A(1,2) + 2.0*C1*DP_C2(Imatrix)*Yfun(Imatrix) - 2*AA1*DP_C2(Imatrix);
        A(1,3) = A(1,3) - 2.0*C1 *YPfun(Imatrix)*Yfun(Imatrix)/C2 + 2*AA1*YPfun(Imatrix)/C2;
        A(1,4) = A(1,4) + 2.0*Yfun(Imatrix);
        A(1,5) = A(1,5) + 2.0*AA1*Yfun(Imatrix);

        %% Ligne 2
        A(2,2) = A(2,2) + 2*C1*C1*DP_C2(Imatrix)*DP_C2(Imatrix) - 2*C1*AA1*DS_C2(Imatrix);
        A(2,3) = A(2,3) - 2*C1/(C2*C2)*AA1*YPfun(Imatrix) - 2*C1.^2/C2*DP_C2(Imatrix)*YPfun(Imatrix) - 2*C1/C2.^3*AA1*YSfun(Imatrix)*CO1;
        A(2,4) = A(2,4) + 2*C1*DP_C2(Imatrix);
        A(2,5) = A(2,5) + 2*C1*AA1*DP_C2(Imatrix);

        %% Ligne 3
        A(3,3) = A(3,3) + 2.0*C1*C1*YPfun(Imatrix)*YPfun(Imatrix)/(C2*C2) - 2*C1*AA1*YSfun(Imatrix)/(C2*C2);
        A(3,4) = A(3,4) - 2*C1*YPfun(Imatrix)/C2;
        A(3,5) = A(3,5) - 2*C1*AA1*YPfun(Imatrix)/C2;

        %% Ligne 4
        A(4,5) = A(4,5) + 2*AA1;
    end
    A(4,4) = 2*(Nfin-Ndep+1);

    %% Copie symétrie
    for i=1:3
        for j=i+1:4
            A(j,i) = A(i,j);
        end
    end

    %% Pivot de gauss
    perm(1:4) = 1:4;
    for K = 1:3
        ma_x = A(K,K);
        lin = K;
        col = K;
        for L = K:4
            for M = K:4
                if abs(A(L,M)) > abs(ma_x)
                    ma_x = A(L,M);
                    col = M;
                    lin = L;
                end
            end
        end
    
        if ~(lin == K)
            for M = 1:5
                R = A(lin,M);
                A(lin,M) = A(K,M);
                A(K,M) = R;
            end
        end
    
        if ~(col == K)
            for M = 1:4
                R = A(M,col);
                A(M,col) = A(M,K);
                A(M,K) = R;
            end
            NN = perm(col);
            perm(col) = perm(K);
            perm(K)= NN;
        end
    
        for L = 1+K:4
            ma_x = A(L,K);
            A(L,K) = 0.0;
            for M = 1+K:5
                A(L,M) = A(L,M) - A(K,M)*ma_x/A(K,K);
            end
        end
    end


    %% Calcul solution-----------------------------------------
    dp(4) = A(4,5)/A(4,4);
    dp(3) = (A(3,5) - dp(4)*A(3,4))/A(3,3);
    dp(2) = (A(2,5)-dp(4)*A(2,4) - dp(3)*A(2,3))/A(2,2);
    dp(1) = A(1,5) - dp(4)*A(1,4) - dp(3)*A(1,3) - dp(2)*A(1,2);
    dp(1) = dp(1)/A(1,1);

    %% remet initial ---------------------------------------
    for L=1:4
        if ~(perm(L)==L)
            NN = perm(L);
            perm(L) = perm(NN);
            perm(NN) = NN;
            Q = dp(L);
            dp(L) = dp(NN);
            dp(NN) = Q;
        end
    end

    %% MAJ les C

    C1 = dp(1) + C1;
    C2 = dp(2) + C2;
    C3 = dp(3) + C3;
    C4 = dp(4) + C4;

    %% Parametres
    gamma = 2.0*Dexp*C2/0.1
    shift = -(Position0-(C2*Point0+C3))*Dexp
end

YF(Ndep:Nfin) = C1*Yfun(Ndep:Nfin)+C4; % 
figure;
plot(Ndep:Nfin,YF(Ndep:Nfin),'r',Ndep:Nfin,Yexp(Ndep:Nfin),'b');

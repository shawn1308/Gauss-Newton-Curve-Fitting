function[fun,funx] = GN_Function_P0(C1,C2,C3,C4)
    global Ythe; global Xthe; global Yexp; global Xexp; 
    global N_DEP; global N_FIN; global Point0;
    funx = zeros(length(N_DEP:N_FIN),1);
    fun = zeros(length(N_DEP:N_FIN),1);
    for i_fn = N_DEP:N_FIN
        xtil = (Xexp(i_fn)-C3)/C2; % Calcul coef Shift + élargissement
        n = max(1,min(floor(xtil)+Point0,length(Ythe)-1)); % Prend entier de cette coef
        d = xtil-n+Point0; % difference entre coef et entier
        idx = i_fn-N_DEP+1;
        fun(idx)= C1*(Ythe(n)*(1-d)+Ythe(n+1)*d)+C4; % interpolation sur la courbe theorique
        funx(idx) = Xexp(i_fn);
    end
end
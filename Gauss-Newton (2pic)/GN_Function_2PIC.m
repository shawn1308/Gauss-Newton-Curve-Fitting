function[fun,funx] = GN_Function_2PIC(C1,C2,C3,C4,C5,C6)
    global Ythe; global Xthe; global Yexp; global Xexp; 
    global N_DEP; global N_FIN; global Point0;
    funx = zeros(length(N_DEP:N_FIN),1);
    fun = zeros(length(N_DEP:N_FIN),1);
    for i_fn = N_DEP:N_FIN
        xtil = (Xexp(i_fn)-C3)/C2;
        n = max(1,min(floor(xtil)+Point0,length(Ythe)-1));
        d = xtil-n+Point0;

        xtil2 = (Xexp(i_fn)-C6)/C2;
        n2 = max(1,min(floor(xtil2)+Point0,length(Ythe)-1));
        d2 = xtil2-n2+Point0;

        idx = i_fn-N_DEP+1;
        fun(idx)= (C1*(Ythe(n)*(1-d)+Ythe(n+1)*d)+C4) + (C5*(Ythe(n2)*(1-d2)+Ythe(n2+1)*d2)+C4);
        funx(idx) = Xexp(i_fn);
    end
end
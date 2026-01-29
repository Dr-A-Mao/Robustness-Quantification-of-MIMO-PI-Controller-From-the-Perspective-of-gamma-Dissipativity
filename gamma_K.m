function gamma_max = gamma_K(A_0,L_K_Omega,S_K_Omega)
    
    P_K = lyap(A_0,eye(size(A_0,1)));

    F_K = [ zeros(size(P_K)) , 2 * (1 - S_K_Omega/L_K_Omega) * P_K;...
            2 * (1 - S_K_Omega/L_K_Omega) * P_K', zeros(size(P_K))];

    gamma_max = max(eig(F_K));

end
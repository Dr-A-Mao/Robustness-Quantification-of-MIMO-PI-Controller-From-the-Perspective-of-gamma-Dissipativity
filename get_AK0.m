function [A_K0,max_real_lambda] = get_AK0(K,J_F0,J_G0)
    % K:(K_P,K_I)的列向量形式
    [num_e,num_u] = size(J_G0);
    K_all = reshape(K,[num_u,2*num_e]);
    K_P = K_all(:,1:num_e);
    K_I = K_all(:,num_e+1:end);
    A_K0 = [J_F0 + J_G0 * K_P,...
                        J_G0 * K_I;...
                        eye(num_e),...
                        zeros(num_e)];
    max_real_lambda = max(real(eig(A_K0)));
end
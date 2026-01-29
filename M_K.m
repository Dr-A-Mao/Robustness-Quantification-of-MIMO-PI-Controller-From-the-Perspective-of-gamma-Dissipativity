function M = M_K(A)
    % 计算特征值
    alphas = real(eig(A)); % 特征值实部
    % 检查是否 Hurwitz 稳定
    if max(alphas) >= 0
        warning(['矩阵不是 Hurwitz 稳定的，' ...
                 'epsilon 将为 0 或负数，M 可能非常大']);
        M = 1e6;
        return;
    end
    % 计算 epsilon
    epsilon = max(alphas);
    % 进行jordan分解
    [P,J] = jordan(A);
    % 计算M
    M = norm(P) * norm(inv(P));
end
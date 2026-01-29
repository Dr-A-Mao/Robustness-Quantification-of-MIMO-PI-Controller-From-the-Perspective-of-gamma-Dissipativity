%% 用于进行Duffingmodel的仿真
clc,clear;close all;
%% 配置Matlab环境参数
set(groot,'defaultLineLineWidth',1);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultAxesFontSize',10);
set(groot,'defaultAxesLabelFontSizeMultiplier',1);
set(groot,'defaultFigurePosition',[600 500 400 300]);
%% 下面开始设置雅可比矩阵
g = 9.81; 
V = 25; 
T = 40; % 40s
chi_c = 0;
gamma_c = pi/12;
chi0 = pi/3; gamma0 = pi/4;
w_chi = 0.15; w_gamma = 0.15;
L_d_chi = 0.1; L_d_gamma = 0.1;
J_F0 = [0,0;0,g/V * sin(gamma_c)];
J_G0 = diag([-g/V,-g/V]);
%% 根据设置的雅可比矩阵求解合适的K_opt
type = 2; % type = 1表示含有共轭项的,type = 2 表示无共轭项的
[num_e,num_u] = size(J_G0);
% !!! K0的初始值是0吗？
zeros_type = 4;
if zeros_type == 0
    K0 = zeros(2*num_e * num_u,1);
elseif zeros_type == 1
    K_P = [2.6347, -5.1510;...
		   -0.5854, 3.1952];
    K_I = [2.4254, -0.7911;...
		   -0.0587, 3.3742];
    K0 = [K_P,K_I];
    K0 = K0(:);
elseif zeros_type == 2
    K_P = [2.5882,0.6054;...
          -0.6054,2.8470];
    K_I = [4.8587,0.1335;...
          -0.1336,4.8588];
    K0 = [K_P,K_I];
    K0 = K0(:);
elseif zeros_type == 3
    K_P = [1.7017,0.5960;...
          -0.5960, 1.9605];
    K_I = [3.5129,0.1857;...
          -0.1857,3.5129];
    K0 = [K_P,K_I];
    K0 = K0(:);
elseif zeros_type == 4
    K_P = [1.6962, 0.5900;...
          -0.5900, 1.9550];
    K_I = [3.4836, 0.1774;...
          -0.1774, 3.4836];
    K0 = [K_P,K_I];
    K0 = K0(:);
end
%% 优化类型
optimization_type = 'fminsearch';
if strcmp(optimization_type,'fminsearch')
    num_MaxIter = 3e4;
    options = optimset('PlotFcns',@optimplotfval,...
                       'MaxIter',num_MaxIter,...
                       'MaxFunEvals',num_MaxIter);
    [K_opt,lamda_min_opt] = fminsearch(@(K)EVP_obj(K,J_F0,J_G0,type),K0,options);
elseif strcmp(optimization_type,'ga')
    options = optimoptions('ga','PlotFcn',@gaplotbestf,...
                           'MaxGenerations',2000,'FunctionTolerance',1e-8);
    [K_opt,lamda_min_opt] = ga(@(K)EVP_obj(K,J_F0,J_G0,type),length(K0),options);
end
K_opt = reshape(K_opt,[num_u,2*num_e]);
K_P_opt = K_opt(:,1:num_e);
K_I_opt = K_opt(:,num_e+1:end);
K_P_opt_new = K_P_opt;
K_I_opt_new = K_I_opt;
%% 仿真结果
out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
disp('-----Simulation completes!-----');
%% 测试A_K(0)的特征值
A_K0 = [J_F0 + J_G0 * K_P_opt_new, J_G0 * K_I_opt_new;eye(num_e),zeros(num_e)];
QK_star = get_QK_star(A_K0);
RK = get_RK(QK_star);
IK = get_IK(QK_star,RK,A_K0);
% 打印其最大特征值的实数部分
disp(['The eigvalues of A_K(0):']);
disp(eig(A_K0));
disp(['R_K:',num2str(RK),',I_K:',num2str(IK)]);
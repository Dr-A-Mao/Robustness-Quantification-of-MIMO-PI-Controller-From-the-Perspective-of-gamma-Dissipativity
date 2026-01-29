%% 用于进行Duffingmodel的仿真
clc,clear;close all;rng(0);
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

%% 绘制仿真结果图Exp1
% 画极点图
figure; hold on;box on; 
real_eig_J = real(eig(A_K0));
imag_eig_J = imag(eig(A_K0));
plot(real_eig_J,imag_eig_J,'r+','LineWidth',1.5); 
xlabel('Real');ylabel('Imag');
max_real_eig_J = max(abs(real_eig_J));
max_imag_eig_J = max(abs(imag_eig_J));
plot([-1.15 * max_real_eig_J, 1.15 * max_real_eig_J],[0,0],...
      'LineWidth',1,'Color','b','LineStyle','--');
plot([0,0],[-1.15 * max_imag_eig_J, 1.15 * max_imag_eig_J],...
      'LineWidth',1,'Color','b','LineStyle','--');
% 画误差收敛图
figure;
hold on; box on;
plot(out.tout, out.e_chi,'b-','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('The error of \chi (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.e_gamma,'b-','LineWidth',1.5,'DisplayName','e_{\gamma}');
plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('The error of \gamma (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.e_chi_dot,'b-','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('The differential error of \chi (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.e_gamma_dot,'b-','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('The differential error of \gamma (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.d_chi,'r--','LineWidth',1.5,'DisplayName','d_{\chi}');
plot(out.tout, out.d_gamma,'b-.','LineWidth',1.5,'DisplayName','d_{\gamma}');
xlabel('t(s)'); ylabel('Disturbances'); legend();

figure;
hold on; box on;
plot(out.tout, out.e_gamma_dot,'b-','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('The differential error of \gamma (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.phi,'b-','LineWidth',1.5,'DisplayName','\phi');
%plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('\phi (rad)'); legend();

figure;
hold on; box on;
plot(out.tout, out.nz,'b-','LineWidth',1.5,'DisplayName','n_z');
%plot(out.tout, zeros(size(out.tout)),'r--','LineWidth',1.5,'DisplayName','origin');
xlabel('t(s)'); ylabel('n_z'); legend();

%% 在K_P_opt和K_I_opt附近施加扰动
varepsilon_array = [-4,-2,-1,0,0.5,0.8,1];
states_array = cell(size(varepsilon_array));
str_array = cell(size(varepsilon_array));
RK_array = zeros(size(varepsilon_array));
IK_array = zeros(size(varepsilon_array));
for k = 1:length(varepsilon_array)
    varepsilon_P = varepsilon_array(k);
    varepsilon_I = varepsilon_array(k);
    % 更新K_P_opt_new, K_I_opt_new
    K_P_opt_new = K_P_opt - varepsilon_P * eye(num_e);
    K_I_opt_new = K_I_opt - varepsilon_I * eye(num_e);
    % 仿真
    out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
    disp(['-----Simulation of ',num2str(varepsilon_P) ,' completes!-----']);
    pause(0.1);
    % 计算对应的RK，IK
    A_K0 = [J_F0 + J_G0 * K_P_opt_new, J_G0 * K_I_opt_new;eye(num_e),zeros(num_e)];
    QK_star = get_QK_star(A_K0);
    RK = get_RK(QK_star);
    IK = get_IK(QK_star,RK,A_K0);
    RK_array(k) = RK;
    IK_array(k) = IK;
    % 计算对应的曲线
    states = zeros(length(out.tout),5);
    % 状态填充
    states(:,1) = out.chi;
    states(:,2) = out.gamma;
    states(:,3) = -out.e_chi_dot;
    states(:,4) = -out.e_gamma_dot;
    states(:,5) = out.tout;
    states_array{k} = states;
end
%% 画4个曲线的对比图
lines_array = {'-', '--', '-.', ':', '-.', '--', ':'};
colors_array = { 
    [0.0000    0.4470    0.7410],...  % 紫色
    [0.8500    0.3250    0.0980],...  % 橙色
    [0.9290    0.6940    0.1250],...  % 黄色
    [0.4940    0.1840    0.5560],...  % 蓝色
    [0.4660    0.6740    0.1880],...  % 绿色
    [0.3010    0.7450    0.9330],...  % 青色
    [0.6350    0.0780    0.1840]  % 红色
};
figure(101);
hold on; box on;
xlabel('t(s)'); ylabel('\chi (rad)');
figure(102);
hold on; box on;
xlabel('t(s)'); ylabel('\gamma (rad)');
figure(103);
hold on; box on;
xlabel('t(s)'); ylabel('The differential of \chi (rad/s)');
figure(104);
hold on; box on;
xlabel('t(s)'); ylabel('The differential of \gamma (rad/s)');
for k = 1:length(varepsilon_array)
    states = states_array{k};
    state_tout = states(:,5);
    state_chi = states(:,1);
    state_gamma = states(:,2);
    state_chi_dot = states(:,3);
    state_gamma_dot = states(:,4);
    % 第1次画图时画参考轨迹
    if k == 1
        figure(101);
        plot(state_tout,chi_c * ones(size(state_tout)),'r--',...
                                 'LineWidth',1.5,'DisplayName','Reference');
        figure(102);
        plot(state_tout,gamma_c * ones(size(state_tout)),'r--',...
                                 'LineWidth',1.5,'DisplayName','Reference');
    end
    % 下面开始画图
    str_array{k} = ['R_K=',num2str(RK_array(k)),...
                    ',I_K=',num2str(IK_array(k))];
    figure(101);
    plot(state_tout,state_chi,'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
    figure(102);
    plot(state_tout,state_gamma,'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
    figure(103);
    plot(state_tout,state_chi_dot,'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
    figure(104);
    plot(state_tout,state_gamma_dot,'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
end
figure(101); legend();
figure(102); legend();
figure(103); legend();
figure(104); legend();
%% 画相应的指标对比图
ITAE_array = zeros(length(varepsilon_array),4);
MO_array = zeros(length(varepsilon_array),4);
PT_array = zeros(length(varepsilon_array),4);
MS_array = zeros(length(varepsilon_array),4);
ST_array = zeros(length(varepsilon_array),4);
for k = 1:length(varepsilon_array)
    states = states_array{k};
    % 第1列是e_chi
    [ITAE_e_chi,MO_e_chi,PT_e_chi,MS_e_chi,ST_e_chi] = cal_key_indictors(states(:,5),chi_c - states(:,1));
    % 第2列是e_gamma
    [ITAE_e_gamma,MO_e_gamma,PT_e_gamma,MS_e_gamma,ST_e_gamma] = cal_key_indictors(states(:,5),gamma_c - states(:,2));
    % 第3列是e_chi dot
    [ITAE_e_chi_dot,MO_e_chi_dot,PT_e_chi_dot,MS_e_chi_dot,ST_e_chi_dot] = cal_key_indictors(states(:,5),- states(:,3));
    % 第4列是e_gamma dot
    [ITAE_e_gamma_dot,MO_e_gamma_dot,PT_e_gamma_dot,MS_e_gamma_dot,ST_e_gamma_dot] = cal_key_indictors(states(:,5),- states(:,4));
    % 赋值给矩阵
    ITAE_array(k,:) = [ITAE_e_chi,ITAE_e_gamma,ITAE_e_chi_dot,ITAE_e_gamma_dot];
    MO_array(k,:) = [MO_e_chi,MO_e_gamma,MO_e_chi_dot,MO_e_gamma_dot];
    PT_array(k,:) = [PT_e_chi,PT_e_gamma,PT_e_chi_dot,PT_e_gamma_dot];
    MS_array(k,:) = [MS_e_chi,MS_e_gamma,MS_e_chi_dot,MS_e_gamma_dot];
    ST_array(k,:) = [ST_e_chi,ST_e_gamma,ST_e_chi_dot,ST_e_gamma_dot];
end
%% 开始画图ITAE,MO,PT
name_array = {'e_{\chi}','e_{\gamma}','dot e_{\chi}','dot e_{\gamma}'};
title_name_ = reordercats(categorical(name_array),name_array);

figure; hold on; box on;% 画ITAE对比
bar(title_name_,ITAE_array,'grouped');
ylabel('ITAE'); legend(str_array);

figure; hold on; box on;% 画MO对比
bar(title_name_,MO_array,'grouped');
ylabel('MO'); legend(str_array);

figure; hold on; box on;% 画PT对比
bar(title_name_,PT_array,'grouped');
ylabel('PT'); legend(str_array);

%% 开始画图ITAE,MO,PT
name_array = {'s_{\chi}','s_{\gamma}'};
title_name_ = reordercats(categorical(name_array),name_array);
% 画MS对比
figure; hold on; box on;
b_d = bar(title_name_,[sqrt(MS_array(:,1).^2 + MS_array(:,3).^2),...
                       sqrt(MS_array(:,2).^2 + MS_array(:,4).^2)],...
          'BarLayout','grouped',...
          'LineWidth',1);
ylabel('MS'); 
legend(str_array);
% 画ST对比
figure; hold on; box on;
b_d = bar(title_name_,[sqrt(ST_array(:,1).^2 + ST_array(:,3).^2),...
                       sqrt(ST_array(:,2).^2 + ST_array(:,4).^2)],...
          'BarLayout','grouped',...
          'LineWidth',1);
ylabel('ST'); 
legend(str_array);
%% R_K - MS, ST
figure; hold on; box on;
plot(RK_array,MS_array(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(RK_array,MS_array(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('R_K'); ylabel('The average value of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(RK_array,MS_array(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(RK_array,MS_array(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('R_K'); ylabel('The average value of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

figure; hold on; box on;
plot(RK_array,ST_array(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(RK_array,ST_array(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('R_K'); ylabel('The standard deviation of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(RK_array,ST_array(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(RK_array,ST_array(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('R_K'); ylabel('The standard deviation of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

%% I_K - MS, ST
MS_array_ = MS_array; ST_array_ = ST_array;
[IK_array,IK_sort_index] = sort(IK_array);
for k = 1:4
    MS_array_(:,k) = MS_array(IK_sort_index,k);
    ST_array_(:,k) = ST_array(IK_sort_index,k);
end

figure; hold on; box on;
plot(IK_array,MS_array_(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(IK_array,MS_array_(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('I_K'); ylabel('The average value of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(IK_array,MS_array_(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(IK_array,MS_array_(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('I_K'); ylabel('The average value of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

figure; hold on; box on;
plot(IK_array,ST_array_(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(IK_array,ST_array_(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('I_K'); ylabel('The standard deviation of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(IK_array,ST_array_(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(IK_array,ST_array_(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('I_K'); ylabel('The standard deviation of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

%% 再对比不同的L_d_chi，w_chi时的曲线
K_I_opt_new = K_I_opt;
K_P_opt_new = K_P_opt;
w_chi_base = w_chi; w_gamma_base = w_gamma;
L_d_chi_base = L_d_chi; L_d_gamma_base = L_d_gamma;
Amplifier_L_array = [0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3]; 
Amplifier_w_array = [0.1,0.15,0.2,0.1,0.15,0.2,0.1,0.15,0.2];

ITAE_array = zeros(size(Amplifier_L_array));
MS_array = zeros(size(Amplifier_L_array));
ST_array = zeros(size(Amplifier_L_array));

tout_array = cell(1,length(Amplifier_L_array));
e_chi_array = cell(1,length(Amplifier_L_array));
e_gamma_array = cell(1,length(Amplifier_L_array));
e_chi_dot_array = cell(1,length(Amplifier_L_array));
e_gamma_dot_array = cell(1,length(Amplifier_L_array));
s_array = cell(1,length(Amplifier_L_array));
s_dot_array = cell(1,length(Amplifier_L_array));
L_d_array = cell(1,length(Amplifier_L_array));
w_array = cell(1,length(Amplifier_L_array));

for k = 1:length(Amplifier_L_array)
    L_d_chi = Amplifier_L_array(k);
    L_d_gamma = Amplifier_L_array(k);
    w_chi = Amplifier_w_array(k);
    w_gamma = Amplifier_w_array(k);
    out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
    % 计算对应的数据
    tout_array{k} = out.tout;
    e_chi_array{k} = out.e_chi;
    e_gamma_array{k} = out.e_gamma;
    e_chi_dot_array{k} = out.e_chi_dot;
    e_gamma_dot_array{k} = out.e_gamma_dot;
    s_array{k} = sqrt((out.e_chi).^2 + (out.e_gamma).^2);
    s_dot_array{k} = sqrt((out.e_chi_dot_array).^2 + (out.e_gamma_dot_array).^2);
    L_d_array{k} = L_d_chi;
    w_array{k} = w_chi;
    % 对应的指标
    [ITAE_s,MO_s,PT_s,MS_s,ST_s] = cal_key_indictors(tout_array{k},s_array{k});
    % 指标
    ITAE_array(k) = ITAE_s;
    MS_array(k) = MS_s;
    ST_array(k) = ST_s;
end
% 画图
figure(201);hold on; box on; xlabel('t (s)'); ylabel('e_{\chi}');
figure(202);hold on; box on; xlabel('t (s)'); ylabel('e_{\gamma}');
figure(203);hold on; box on; xlabel('t (s)'); ylabel('dot e_{\chi}');
figure(204);hold on; box on; xlabel('t (s)'); ylabel('dot e_{\gamma}');
lines_array = {'-', '--', '-.','-', '--','-.','-','--','-.'};
colors_array = { 
    [0.0000    0.4470    0.7410],...  % 紫色
    [0.8500    0.3250    0.0980],...  % 橙色
    [0.9290    0.6940    0.1250],...  % 黄色
    [0.4940    0.1840    0.5560],...  % 蓝色
    [0.4660    0.6740    0.1880],...  % 绿色
    [0.3010    0.7450    0.9330],...  % 青色
    [0.6350    0.0780    0.1840],...  % 红色
    'r','b'};
for k = 1:length(Amplifier_L_array)
    L_d = L_d_array{k};
    w_d = w_array{k};
    DisplayName_str = ['L_d=',num2str(L_d(1)),...
                       ',\omega_d=',num2str(w_d(1))];
    figure(201);
    plot(tout_array{k},e_chi_array{k},...
             'color',colors_array{k},...
             'linewidth',1.5,...
             'linestyle',lines_array{k},...
             'DisplayName',DisplayName_str);
    figure(202);
    plot(tout_array{k},e_gamma_array{k},...
             'color',colors_array{k},...
             'linewidth',1.5,...
             'linestyle',lines_array{k},...
             'DisplayName',DisplayName_str);
    figure(203);
    plot(tout_array{k},e_chi_dot_array{k},...
             'color',colors_array{k},...
             'linewidth',1.5,...
             'linestyle',lines_array{k},...
             'DisplayName',DisplayName_str);
    figure(204);
    plot(tout_array{k},e_gamma_dot_array{k},...
             'color',colors_array{k},...
             'linewidth',1.5,...
             'linestyle',lines_array{k},...
             'DisplayName',DisplayName_str);
end
figure(201);legend();
figure(202);legend();
figure(203);legend();
figure(204);legend();
%% 画ITAE, MO, PT的图
figure; hold on; box on;legend();
xlabel('w (rad/s)'); ylabel('ITAE');
plot(Amplifier_w_array(1:3),ITAE_array(1:3),...
               'color',colors_array{1},'LineStyle',lines_array{1},...
               'Marker','s','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(1))]);
plot(Amplifier_w_array(4:6),ITAE_array(4:6),...
               'color',colors_array{2},'LineStyle',lines_array{2},...
               'Marker','d','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(4))]);
plot(Amplifier_w_array(7:9),ITAE_array(7:9),...
               'color',colors_array{3},'LineStyle',lines_array{3},...
               'Marker','o','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(7))]);

figure; hold on; box on;legend();
xlabel('w (rad/s)'); ylabel('MS');
plot(Amplifier_w_array(1:3),MS_array(1:3),...
               'color',colors_array{1},'LineStyle',lines_array{1},...
               'Marker','s','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(1))]);
plot(Amplifier_w_array(4:6),MS_array(4:6),...
               'color',colors_array{2},'LineStyle',lines_array{2},...
               'Marker','d','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(4))]);
plot(Amplifier_w_array(7:9),MS_array(7:9),...
               'color',colors_array{3},'LineStyle',lines_array{3},...
               'Marker','o','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(7))]);

figure; hold on; box on;legend();
xlabel('w (rad/s)'); ylabel('ST');
plot(Amplifier_w_array(1:3),ST_array(1:3),...
               'color',colors_array{1},'LineStyle',lines_array{1},...
               'Marker','s','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(1))]);
plot(Amplifier_w_array(4:6),ST_array(4:6),...
               'color',colors_array{2},'LineStyle',lines_array{2},...
               'Marker','d','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(4))]);
plot(Amplifier_w_array(7:9),ST_array(7:9),...
               'color',colors_array{3},'LineStyle',lines_array{3},...
               'Marker','o','Linewidth',1.5,...
               'DisplayName',['L_d=',num2str(Amplifier_L_array(7))]);




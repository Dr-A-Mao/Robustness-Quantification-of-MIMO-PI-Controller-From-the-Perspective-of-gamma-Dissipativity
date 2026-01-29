%% 用于进行Duffingmodel的仿真
clc,clear;close all;rng(0);
%% 配置Matlab环境参数
set(groot,'defaultLineLineWidth',1);
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultAxesFontSize',10);
set(groot,'defaultAxesLabelFontSizeMultiplier',1);
set(groot,'defaultFigurePosition',[600 500 400 300]);
%% 基本参数
markers_array = {'o','d','s','*','+'};
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
%% 下面开始设置雅可比矩阵
g = 9.81; 
V = 25; 
T = 40; % 40s
% 初始条件和参考量
chi_c = 0; gamma_c = pi/12;
chi0 = pi/3; gamma0 = pi/4;
% 扰动项
w_chi = 0.15; w_gamma = 0.15;
L_d_chi = 0.1; L_d_gamma = 0.1;
% 控制器参数
K_P0 = [1.6962, 0.5900;...
        -0.5900, 1.9550];
K_I0 = [3.4836, 0.1774;...
        -0.1774, 3.4836];

num_e = size(K_P0,1);
varepsilon_array = [-4,-2,-1,0,0.5,0.8,1];

J_F0 = [0,0;0,g/V * sin(gamma_c)];
J_G0 = diag([-g/V,-g/V]);

% e = [e_chi,e_gamma]', u = [phi,nz]'
J_F = @(e,u)[0,0; 0, g/V * sin(gamma_c - e(2))];
J_G = @(e,u)[-g/V /(cos(u(1)))^2, 0; g * u(2)/V * sin(u(1)), -g/V * cos(u(1)) ];

A_K = @(K_P, K_I,e,u) [J_F(e,u)  + J_G(e,u) * K_P, J_G(e,u) * K_I; ...
                       eye(length(e)), zeros(length(e))];
A_0 = @(K_P, K_I)[J_F0 + J_G0 * K_P, J_G0 * K_I;...
                  eye(size(K_P,2)), zeros(size(K_P,2))];
L_K = @(K_P, K_I,e,u) [max(real(eig(A_K(K_P, K_I,e,u)))) +...
                       norm(A_K(K_P, K_I,e,u) - A_0(K_P, K_I)) * (M_K(A_0(K_P, K_I)))^2];
S_K = @(K_P, K_I,e,u)[norm(A_K(K_P, K_I,e,u) - A_0(K_P, K_I)) * (M_K(A_0(K_P, K_I)))^2];
%% 绘制在K=K0时的e_chi, e_gamma的对应的L_K
if false
    disp('------算法最优参数搜索------');
    for k = 1:length(varepsilon_array)
        varepsilon_P = varepsilon_array(k);
        varepsilon_I = varepsilon_array(k);
        disp(['正在进行varepsilon_P=',num2str(varepsilon_P),'时的测试！']);
        if varepsilon_P == 0
            continue;
        end
        % 更新K_P_opt_new, K_I_opt_new
        K_P_opt_new = K_P0 - varepsilon_P * eye(num_e);
        K_I_opt_new = K_I0 - varepsilon_I * eye(num_e);
        % 仿真计算对应的应力图
        e_chi_array = -pi/3:0.05:pi/3;
        e_gamma_array = -pi/6:0.01:pi/6;
        num_e_chi = length(e_chi_array);
        num_e_gamma = length(e_gamma_array);
        [E_CHI,E_GAMMA] = meshgrid(e_chi_array,e_gamma_array);
        L_E = zeros(num_e_gamma,num_e_chi);
        S_E = zeros(num_e_gamma,num_e_chi);
        for i = 1:num_e_chi
            for j = 1:num_e_gamma
                e_chi = e_chi_array(i);
                e_gamma = e_gamma_array(j);
                clc;disp(['------已完成:',num2str(100 * (num_e_gamma * (i-1) + j)/(num_e_chi * num_e_gamma)),'%------']);
                L_E(j,i) = L_K(K_P_opt_new,K_I_opt_new,[e_chi,e_gamma]',[0,1]');
                S_E(j,i) = L_E(j,i) - max(real(eig(A_K(K_P_opt_new,K_I_opt_new,[e_chi,e_gamma]',[0,1]'))));
            end
        end
        save(['L_K_',num2str(k),'.mat']);
    end
end
%% 跑一个小实验
if false
    try
        load('L_K.mat');
    catch
        e_chi_array = -pi/3:0.05:pi/3;
        e_gamma_array = -pi/6:0.01:pi/6;
        num_e_chi = length(e_chi_array);
        num_e_gamma = length(e_gamma_array);
        [E_CHI,E_GAMMA] = meshgrid(e_chi_array,e_gamma_array);
        L_E = zeros(num_e_gamma,num_e_chi);
        for i = 1:num_e_chi
            for j = 1:num_e_gamma
                e_chi = e_chi_array(i);
                e_gamma = e_gamma_array(j);
                clc;disp(['------已完成:',num2str(100 * (num_e_gamma * (i-1) + j)/(num_e_chi * num_e_gamma)),'%------']);
                L_E(j,i) = L_K(K_P0,K_I0,[e_chi,e_gamma]',[0,1]');
            end
        end
        % 存储数据
        save('L_K.mat');
    end
    
    S_E = zeros(num_e_gamma,num_e_chi);
    for i = 1:num_e_chi
        for j = 1:num_e_gamma
            e_chi = e_chi_array(i);
            e_gamma = e_gamma_array(j);
            clc;disp(['------已完成:',num2str(100 * (num_e_gamma * (i-1) + j)/(num_e_chi * num_e_gamma)),'%------']);
            S_E(j,i) = L_E(j,i) - max(real(eig(A_K(K_P0, K_I0,[e_chi,e_gamma]',[0,1]'))));
        end
    end

end

%% 绘制gamma_K_star曲线图，展示不同参数下的耗散域
num_L_K_conditions = length(varepsilon_array); 
vmin = -1; vmax = 0.2;
resolution_array = [0.01, 0.005, 0.004, 0.003, 0.003, 0.003, 0.003];
gamma_K_star_matrix = zeros(num_L_K_conditions,4);
L_omega_array = zeros(num_L_K_conditions,1);
fig_handles = cell(1,num_L_K_conditions);
for L_K_condition = 1:num_L_K_conditions
    % 创建新图
    fig_handles{L_K_condition} = figure;
    % 导入数据
    load(['L_K_',num2str(L_K_condition),'.mat']);
    % 设计插值函数
    L = @(e_chi,e_gamma)interp2(E_CHI,E_GAMMA,L_E,e_chi,e_gamma,'cubic');
    S = @(e_chi,e_gamma)interp2(E_CHI,E_GAMMA,S_E,e_chi,e_gamma,'cubic');
    % 画一个大图
    hold on; grid on; box on;
    h_mesh = mesh(E_CHI,E_GAMMA,L_E);
    colorbar; 
    shading interp;  % 平滑着色（可选）
    caxis([vmin vmax]);% 核心：设置统一的颜色轴范围
    xlabel('e_{\chi}(rad)');ylabel('e_{\gamma}(rad)');zlabel('L_K');
    xlim([e_chi_array(1),e_chi_array(end)]); ylim([e_gamma_array(1),e_gamma_array(end)]);
    % 画等高线
    L_K_Omega = -0.01; 
    S_K_Omega = 1.5;
    [eh_row,eh_col]= find(abs(L_E - L_K_Omega)<=resolution_array(L_K_condition));
    zeros_line_x = zeros(1,length(eh_row));
    zeros_line_y = zeros(1,length(eh_row));
    zeros_line_z = zeros(1,length(eh_row));
    for k = 1:length(eh_row)
        LE_i = eh_row(k);
        LE_j = eh_col(k);
        
        zeros_line_x(k) = e_chi_array(LE_j);
        zeros_line_y(k) = e_gamma_array(LE_i);
        zeros_line_z(k) = L_E(LE_i,LE_j);
    end
    L_omega_array(L_K_condition) = max(zeros_line_y) - min(zeros_line_y);
    h_zero_line = plot3(zeros_line_x, zeros_line_y, zeros_line_z, ...
               'DisplayName','zero line',...
               'Marker' ,'.','Color','r','LineWidth',1, 'LineStyle','none');
    h_origin = plot3(0,0,0,'Color','k','LineStyle' ,'none',...
                          'LineWidth',1.5,'Marker','o','DisplayName','Origin');
    legend([h_zero_line,h_origin]);
    % 在当前基础上绘制区域框线
    num_bound_conditions = 4;
    for bound_condition = 1:num_bound_conditions
        % 获得边界条件
        if bound_condition == 1
            e_chi_bound = [-0.7, 0.7];
            e_gamma_bound = [-0.3, 0.3];
            bound_color = 'm';
            bound_style = '-.';
        elseif bound_condition == 2
            e_chi_bound = [-0.5, 0.5];
            e_gamma_bound = [-0.22, 0.22];
            bound_color = 'b';
            bound_style = '-.';
        elseif bound_condition == 3
            e_chi_bound = [-0.25, 0.25];
            e_gamma_bound = [-0.15, 0.15];
            bound_color = 'g';
            bound_style = '-.';
        elseif bound_condition == 4
            e_chi_bound = [-0.1, 0.1];
            e_gamma_bound = [-0.06, 0.06];
            bound_color = 'r';
            bound_style = '-.';    
        end
        % 计算在某个区域\Omega内的能量增益\gamma_K_star
        num_area_e_chi = 40;
        num_area_e_gamma = 20;
        L_matrix = zeros(num_area_e_gamma,num_area_e_chi);
        S_matrix = zeros(num_area_e_gamma,num_area_e_chi);
        e_chi_bound_array = linspace(e_chi_bound(1),e_chi_bound(2),num_area_e_chi);
        e_gamma_bound_array = linspace(e_gamma_bound(1),e_gamma_bound(2),num_area_e_gamma);
        for i = 1:num_area_e_chi
            for j = 1:num_area_e_gamma
                e_chi = e_chi_bound_array(i);
                e_gamma = e_gamma_bound_array(j);
                L_matrix(j,i) = L(e_chi,e_gamma);
                S_matrix(j,i) = S(e_chi,e_gamma);
            end
        end
        L_K_Omega = max(max(L_matrix));
        S_K_Omega = max(max(S_matrix));
        % 计算矩形区域内的gamma_max
        gamma_max = gamma_K(A_0(K_P_opt_new, K_I_opt_new),L_K_Omega,S_K_Omega);
        % 赋值
        gamma_K_star_matrix(L_K_condition, bound_condition) = gamma_max;
        % 画方框图
        plot3([e_chi_bound(1),e_chi_bound(2),e_chi_bound(2),e_chi_bound(1),e_chi_bound(1)],...
              [e_gamma_bound(2),e_gamma_bound(2),e_gamma_bound(1),e_gamma_bound(1),e_gamma_bound(2)],...
              [0,0,0,0,0],'Color',bound_color, ...
              'LineStyle',bound_style,'LineWidth',1.5, ...
              'DisplayName',['\gamma_{K}(\Omega_',num2str(bound_condition),')=',num2str(gamma_max,'%.2f')]);
    end
    % 画完框线后开始画轨迹
    if false % 图太难看了
        % 画各个line
        num_conditions = 4;
        handles_array = [];
        for initial_condition = 1:num_conditions
            % 声明图的绘制形式
            if initial_condition == 1
                chi0 = pi/3; gamma0 = pi/4;
                linestyle = '--'; linecolor = 'b';
            elseif initial_condition == 2
                chi0 = pi/4; gamma0 = pi/4;
                linestyle = '-'; linecolor = 'g';
            elseif initial_condition == 3
                chi0 = pi/4; gamma0 = -pi/12;
                linestyle = '-.'; linecolor = 'm';
            elseif initial_condition == 4
                chi0 = pi/3; gamma0 = -pi/12;
                linestyle = '-.'; linecolor = 'c';
            elseif initial_condition == 5
                chi0 = pi/3; gamma0 = pi/12;
                linestyle = '-.'; linecolor = '#EDB120';
            end
            % 进行仿真
            out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
            % 进行相图绘制
            plot3(out.e_chi,out.e_gamma,zeros(size(out.e_chi)), ...
               'LineStyle',linestyle,'Color',linecolor,'Linewidth',1.25,...
               'DisplayName',['\chi(0):',num2str(chi0),',\gamma(0):',num2str(gamma0)])
        end

    end
end
return;

%% 下面配置不同边界下的区域计算对应的gamma_max
% num_bound_conditions = 4;
% for bound_condition = 1:num_bound_conditions
%     % 获得边界条件
%     if bound_condition == 1
%         e_chi_bound = [-0.7, 0.7];
%         e_gamma_bound = [-0.3, 0.3];
%         bound_color = 'm';
%         bound_style = '-.';
%     elseif bound_condition == 2
%         e_chi_bound = [-0.5, 0.5];
%         e_gamma_bound = [-0.22, 0.22];
%         bound_color = 'b';
%         bound_style = '-.';
%     elseif bound_condition == 3
%         e_chi_bound = [-0.25, 0.25];
%         e_gamma_bound = [-0.15, 0.15];
%         bound_color = 'g';
%         bound_style = '-.';
%     elseif bound_condition == 4
%         e_chi_bound = [-0.1, 0.1];
%         e_gamma_bound = [-0.06, 0.06];
%         bound_color = 'r';
%         bound_style = '-.';    
%     end
%     % 计算在某个区域\Omega内的能量增益\gamma_K_star
%     num_area_e_chi = 40;
%     num_area_e_gamma = 20;
%     L_matrix = zeros(num_area_e_gamma,num_area_e_chi);
%     S_matrix = zeros(num_area_e_gamma,num_area_e_chi);
%     e_chi_bound_array = linspace(e_chi_bound(1),e_chi_bound(2),num_area_e_chi);
%     e_gamma_bound_array = linspace(e_gamma_bound(1),e_gamma_bound(2),num_area_e_gamma);
%     for i = 1:num_area_e_chi
%         for j = 1:num_area_e_gamma
%             e_chi = e_chi_bound_array(i);
%             e_gamma = e_gamma_bound_array(j);
%             L_matrix(j,i) = L(e_chi,e_gamma);
%             S_matrix(j,i) = S(e_chi,e_gamma);
%         end
%     end
%     L_K_Omega = max(max(L_matrix));
%     S_K_Omega = max(max(S_matrix));
%     % 计算矩形区域内的gamma_max
%     gamma_max = gamma_K(A_0(K_P0, K_I0),L_K_Omega,S_K_Omega);
%     % 画方框图
%     plot3([e_chi_bound(1),e_chi_bound(2),e_chi_bound(2),e_chi_bound(1),e_chi_bound(1)],...
%           [e_gamma_bound(2),e_gamma_bound(2),e_gamma_bound(1),e_gamma_bound(1),e_gamma_bound(2)],...
%           [0,0,0,0,0],'Color',bound_color, ...
%           'LineStyle',bound_style,'LineWidth',1.5, ...
%           'DisplayName',['\gamma_{K^*}(\Omega)=',num2str(gamma_max,'%.2f')]);
% end

%% 绘制状态轨迹
% out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
% plot3( out.e_chi,out.e_gamma,zeros(size(out.e_chi)), ...
%        'LineStyle','-','Color','k','Linewidth',1,...
%        'DisplayName',['e(t)']);
% 
% disp(['gamma_K_star:',num2str(gamma_max)]);
% disp('-----Simulation completes!-----');
% return;

%% 不同初始状态下采用K^*时的状态轨迹对比
figure;
% 画一个大图
hold on; grid on; box on;
h_mesh = mesh(E_CHI,E_GAMMA,L_E);
colorbar; load('L_K_4.mat');
xlabel('e_{\chi}(rad)');ylabel('e_{\gamma}(rad)');zlabel('L_K');
xlim([e_chi_array(1),e_chi_array(end)]); 
ylim([e_gamma_array(1),e_gamma_array(end)]);
% 画等高线
L_K_Omega = -0.01; 
S_K_Omega = 1.5;
[eh_row,eh_col]= find(abs(L_E - L_K_Omega)<=resolution_array(L_K_condition));
zeros_line_x = zeros(1,length(eh_row));
zeros_line_y = zeros(1,length(eh_row));
zeros_line_z = zeros(1,length(eh_row));
for k = 1:length(eh_row)
    LE_i = eh_row(k);
    LE_j = eh_col(k);
    
    zeros_line_x(k) = e_chi_array(LE_j);
    zeros_line_y(k) = e_gamma_array(LE_i);
    zeros_line_z(k) = L_E(LE_i,LE_j);
end
% 画各个line
num_conditions = 4;
handles_array = [];
for initial_condition = 1:num_conditions
    % 声明图的绘制形式
    if initial_condition == 1
        chi0 = pi/3; gamma0 = pi/4;
        linestyle = '--'; linecolor = 'b';
    elseif initial_condition == 2
        chi0 = pi/4; gamma0 = pi/4;
        linestyle = '-'; linecolor = 'g';
    elseif initial_condition == 3
        chi0 = pi/4; gamma0 = -pi/12;
        linestyle = '-.'; linecolor = 'm';
    elseif initial_condition == 4
        chi0 = pi/3; gamma0 = -pi/12;
        linestyle = '-.'; linecolor = 'c';
    elseif initial_condition == 5
        chi0 = pi/3; gamma0 = pi/12;
        linestyle = '-.'; linecolor = '#EDB120';
    end
    % 进行仿真
    out = sim('sim_model.slx', [0,T]); % Ts内的仿真时间
    % 进行相图绘制
    handles_array =[ handles_array, plot3( ...
        out.e_chi,out.e_gamma,zeros(size(out.e_chi)), ...
       'LineStyle',linestyle,'Color',linecolor,'Linewidth',1.25,...
       'DisplayName',['\chi(0):',num2str(chi0),',\gamma(0):',num2str(gamma0)]) ];
end

handles_array = [handles_array,plot3(zeros_line_x, zeros_line_y, zeros_line_z, ...
           'DisplayName','zero line',...
           'Marker' ,'.','Color','r','LineWidth',1, 'LineStyle','none')];
handles_array = [handles_array,plot3(0,0,0,'Color','k','LineStyle' ,'none',...
                      'LineWidth',1.5,'Marker','o','DisplayName','Origin')];
legend(handles_array);
return;

%% 在K_P_opt和K_I_opt附近施加扰动
states_array = cell(size(varepsilon_array));
str_array = cell(size(varepsilon_array));
RK_array = zeros(size(varepsilon_array));
IK_array = zeros(size(varepsilon_array));
for k = 1:length(varepsilon_array)
    varepsilon_P = varepsilon_array(k);
    varepsilon_I = varepsilon_array(k);
    % 更新K_P_opt_new, K_I_opt_new
    K_P_opt_new = K_P0 - varepsilon_P * eye(num_e);
    K_I_opt_new = K_I0 - varepsilon_I * eye(num_e);
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

%% 计算指标
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
%% 绘制响应曲线
figure(101);
hold on; box on;
xlabel('t(s)'); ylabel(' ||(e_\chi, e_{\gamma})|| (rad)');
figure(102);
hold on; box on;
xlabel('t(s)'); ylabel(' ||(dot e_\chi, dot e_{\gamma})|| (rad/s)');
for k = 1:length(varepsilon_array)
    states = states_array{k};
    state_tout = states(:,5);
    state_chi = states(:,1);
    state_gamma = states(:,2);
    state_chi_dot = states(:,3);
    state_gamma_dot = states(:,4);

    % 下面开始画图
    str_array{k} = ['W_K=',num2str( L_omega_array(k) ),...
                    ', \gamma_K(\Omega_3)=',num2str( gamma_K_star_matrix(k,3)),...
                    ', \gamma_K(\Omega_4)=',num2str( gamma_K_star_matrix(k,4))];
    figure(101);
    plot(state_tout, ...
        sqrt((chi_c - state_chi).^2 + (gamma_c - state_gamma).^2),...
        'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
    figure(102);
    plot(state_tout, ...
        sqrt((- state_chi_dot).^2 + (- state_gamma_dot).^2), ...
        'Color',colors_array{k},...
        'Linewidth',1.5,'DisplayName',str_array{k},...
                                       'LineStyle',lines_array{k});
end
figure(101); legend();
figure(102); legend();
return;

%% 绘制RK和ITAE的关系
figure;
hold on; box on;
xlabel('R_K'); ylabel('ITAE of e_{\chi} and e_{\gamma} (rad)');
plot(RK_array',ITAE_array(:,1), 'Displayname', 'e_{\chi}',...
     'Color',colors_array{1},'LineStyle',lines_array{1}, ...
     'Linewidth',1.5,'Marker',markers_array{1});
plot(RK_array',ITAE_array(:,2), 'Displayname', 'e_{\gamma}',...
     'Color',colors_array{2},'LineStyle',lines_array{2}, ...
     'Linewidth',1.5,'Marker',markers_array{2}); legend();
figure;
hold on; box on;
xlabel('R_K'); ylabel('ITAE of deviations of e_{\chi} and e_{\gamma} (rad)');
plot(RK_array',ITAE_array(:,3), 'Displayname', 'dot e_{\chi}',...
     'Color',colors_array{3},'LineStyle',lines_array{3}, ...
     'Linewidth',1.5,'Marker',markers_array{3});
plot(RK_array',ITAE_array(:,4), 'Displayname', 'dot e_{\gamma}',...
     'Color',colors_array{4},'LineStyle',lines_array{4}, ...
     'Linewidth',1.5,'Marker',markers_array{4});legend();

figure; hold on; box on;
plot(RK_array,ST_array(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(RK_array,ST_array(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('R_K'); ylabel('The standard deviations of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(RK_array,ST_array(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(RK_array,ST_array(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('R_K'); ylabel('The standard deviations of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

return;

%% 绘制gamma_K和ITAE的关系
% figure;
% hold on; box on;
% gamma_K_2_array = gamma_K_star_matrix(:,2);
% xlabel('gamma_K(\Omega_2)'); ylabel('ITAE of e_{\chi} and e_{\gamma} (rad)');
% plot(gamma_K_2_array,ITAE_array(:,1), 'Displayname', 'e_{\chi}',...
%      'Color',colors_array{1},'LineStyle',lines_array{1}, ...
%      'Linewidth',1.5,'Marker',markers_array{1});
% plot(gamma_K_2_array,ITAE_array(:,2), 'Displayname', 'e_{\gamma}',...
%      'Color',colors_array{2},'LineStyle',lines_array{2}, ...
%      'Linewidth',1.5,'Marker',markers_array{2}); legend();
% figure;
% hold on; box on;
% xlabel('gamma_K(\Omega_2)'); ylabel('ITAE of deviations of e_{\chi} and e_{\gamma} (rad)');
% plot(gamma_K_2_array,ITAE_array(:,3), 'Displayname', 'dot e_{\chi}',...
%      'Color',colors_array{3},'LineStyle',lines_array{3}, ...
%      'Linewidth',1.5,'Marker',markers_array{3});
% plot(gamma_K_2_array,ITAE_array(:,4), 'Displayname', 'dot e_{\gamma}',...
%      'Color',colors_array{4},'LineStyle',lines_array{4}, ...
%      'Linewidth',1.5,'Marker',markers_array{4});legend();
% 
% figure; hold on; box on;
% plot(gamma_K_2_array,ST_array(:,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
% plot(gamma_K_2_array,ST_array(:,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
% xlabel('gamma_K(\Omega_2)'); ylabel('The standard deviations of e_{\chi} and e_{\gamma}(rad)');legend();
% 
% figure; hold on; box on;
% plot(gamma_K_2_array,ST_array(:,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
% plot(gamma_K_2_array,ST_array(:,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
% xlabel('gamma_K(\Omega_2)'); ylabel('The standard deviations of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();

%% 绘制width_K和ITAE的关系
figure;
hold on; box on;
[W_K_array,W_K_array_sort_index] = sort(L_omega_array);
xlabel('W_K'); ylabel('ITAE of e_{\chi} and e_{\gamma} (rad)');
plot(W_K_array,ITAE_array(W_K_array_sort_index,1), 'Displayname', 'e_{\chi}',...
     'Color',colors_array{1},'LineStyle',lines_array{1}, ...
     'Linewidth',1.5,'Marker',markers_array{1});
plot(W_K_array,ITAE_array(W_K_array_sort_index,2), 'Displayname', 'e_{\gamma}',...
     'Color',colors_array{2},'LineStyle',lines_array{2}, ...
     'Linewidth',1.5,'Marker',markers_array{2}); legend();
figure;
hold on; box on;
xlabel('W_K'); ylabel('ITAE of deviations of e_{\chi} and e_{\gamma} (rad)');
plot(W_K_array,ITAE_array(W_K_array_sort_index,3), 'Displayname', 'dot e_{\chi}',...
     'Color',colors_array{3},'LineStyle',lines_array{3}, ...
     'Linewidth',1.5,'Marker',markers_array{3});
plot(W_K_array,ITAE_array(W_K_array_sort_index,4), 'Displayname', 'dot e_{\gamma}',...
     'Color',colors_array{4},'LineStyle',lines_array{4}, ...
     'Linewidth',1.5,'Marker',markers_array{4});legend();

figure; hold on; box on;
plot(W_K_array,ST_array(W_K_array_sort_index,1),'Color','#D95319','LineStyle','--','Marker','+','LineWidth',1.5,'DisplayName','e_{\chi}');
plot(W_K_array,ST_array(W_K_array_sort_index,2),'Color','#0072BD','LineStyle','-.','Marker','*','LineWidth',1.5,'DisplayName','e_{\gamma}');
xlabel('W_K'); ylabel('The standard deviations of e_{\chi} and e_{\gamma}(rad)');legend();

figure; hold on; box on;
plot(W_K_array,ST_array(W_K_array_sort_index,3),'Color','#4DBEEE','LineStyle','--','Marker','o','LineWidth',1.5,'DisplayName','dot e_{\chi}');
plot(W_K_array,ST_array(W_K_array_sort_index,4),'Color','#A2142F','LineStyle','-.','Marker','s','LineWidth',1.5,'DisplayName','dot e_{\gamma}');
xlabel('W_K'); ylabel('The standard deviations of dot e_{\chi} and dot e_{\gamma}(rad/s)');legend();


%% 绘制RK和gamma_K和Width_K的关系
figure;
hold on; box on;
xlabel('R_K'); ylabel('Width of L_K(\Omega) <0: W_K');
plot(RK_array(1:length(L_omega_array))',L_omega_array,'r-','LineWidth',1.5,'Marker','+','MarkerEdgeColor','b');

figure;
hold on; box on; 
xlabel('R_K'); ylabel('\gamma_K(\Omega_i) at different areas');
for Omega_index = 1:4
    plot(RK_array(1:length(L_omega_array))',gamma_K_star_matrix(:,Omega_index), 'LineWidth',1.5, ...
         'Marker',markers_array{Omega_index},'Color',colors_array{Omega_index},'LineStyle',lines_array{Omega_index}, ...
         'DisplayName',['\Omega_',num2str(Omega_index)]);
end
legend();

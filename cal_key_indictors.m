function [ITAE,max_value,up_time,stable_value,std_value] = cal_key_indictors(t,e)
% 阈值,当误差稳定小于这个阈值时任务渐进稳定
%epsilon = 0.001;
%epsilon_dot = 0.001;
%e_diff = [0;diff(e)];
t_diff = [eps;eps + diff(t)];
%e_dot = e_diff./t_diff;
% 输入的t和e均为列向量
% ITAE
ITAE = sum(t_diff .* abs(e))/sum(t_diff);
% up_time 
[max_value,max_time_index] = max(e);
up_time = t(max_time_index);
num_e = length(e);
e_new = e(round(0.75 * num_e):end);
stable_value = mean(e_new);
std_value = std(e_new);
% % 若0.9T后的区域都小于epsilon则稳定
% if max(abs(e(find(t>=0.9*t(end))))) <= epsilon
%     stable_time_points = find(abs(e) <= epsilon &...
%                               abs(e_dot) <= epsilon_dot);
%     if isempty(stable_time_points)
%         stable_time = inf;
%         stable_error = inf;
%         disp('需要重新配置参数epsilon_dot!');
%         return;
%     end
%     stable_time = t(stable_time_points(1));
%     stable_error = e(end);
% else
%     stable_time = inf;
%     stable_error = inf;
% end
end
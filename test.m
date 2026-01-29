
% 生成示例数据
[X,Y] = meshgrid(1:5, 1:4);  % 创建网格
Z = [10 15 13 17 14;        % Z值表示每个柱子的高度
     8  12 19 11 16;
     5  9  14 10 7;
     11 13 8  12 9];

% 绘制基本三维柱状图
figure;
bar3(Z);

% 自定义图表外观
colormap jet;               % 设置颜色映射
colorbar;                   % 添加颜色条
shading interp;             % 平滑着色

% 添加标题和坐标轴标签
title('三维柱状图示例');
xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');

% 调整视角和光照
view(30, 40);               % 设置视角
lightangle(45, 30);         % 添加光源
lighting gouraud;           % 设置光照模式

% 为每个柱子添加数值标签
for i = 1:size(Z, 1)
    for j = 1:size(Z, 2)
        text(j, i, Z(i,j)+0.5, num2str(Z(i,j)), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom');
    end
end
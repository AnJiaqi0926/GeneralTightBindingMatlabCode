function PlotGrapheneStructure(a1, a2, a, atom_position, Nx, Ny)
    % 绘制石墨烯结构图函数
    %
    % 输入参数:
    %   a1, a2 - 超胞基矢
    %   a - 晶格常数
    %   atom_position - 原子位置矩阵
    %   Nx, Ny - 超胞尺寸
    
    % 创建新图形窗口
    figure('Name', sprintf('%dx%d 石墨烯超胞结构', Nx, Ny), ...
           'Position', [100, 100, 800, 600]);
    
    % 提取原子位置信息
    x = atom_position(:, 1);          % x坐标
    y = atom_position(:, 2);          % y坐标
    sublattice = atom_position(:, 5); % 子晶格标识 (1=A, 2=B)
    
    % 分离A、B子晶格原子
    idxA = (sublattice == 1);  % A子晶格原子索引
    idxB = (sublattice == 2);  % B子晶格原子索引
    
    % 绘制图形
    hold on;
    
    % ================== 绘制超胞边界 ==================
    % 计算超胞顶点坐标
    vertices = [0, 0; 
                a1(1), a1(2);
                a1(1) + a2(1), a1(2) + a2(2);
                a2(1), a2(2);
                0, 0];
    
    % 绘制超胞边界
    plot(vertices(:, 1), vertices(:, 2), 'k-', 'LineWidth', 2);
    
    % ================== 绘制原子 ==================
    % 绘制A子晶格原子 (红色)
    scatter(x(idxA), y(idxA), 100, 'filled', ...
            'MarkerFaceColor', [0.8, 0.2, 0.2], ...  % 红色
            'MarkerEdgeColor', 'k');
    
    % 绘制B子晶格原子 (蓝色)
    scatter(x(idxB), y(idxB), 100, 'filled', ...
            'MarkerFaceColor', [0.2, 0.2, 0.8], ...  % 蓝色
            'MarkerEdgeColor', 'k');
    
    % ================== 绘制原子键 ==================
    % 仅对小型超胞绘制原子键 (避免性能问题)
    if Nx <= 10 && Ny <= 10
        d_nn = norm(a/sqrt(3));  % 最近邻距离
        
        % 遍历所有原子对
        for i = 1:size(atom_position, 1)
            for j = i+1:size(atom_position, 1)
                % 计算原子间距离
                dx = x(j) - x(i);
                dy = y(j) - y(i);
                dist = sqrt(dx^2 + dy^2);
                
                % 如果是最近邻原子对，绘制键
                if abs(dist - d_nn) < 0.1  % 允许0.1Å的容差
                    plot([x(i), x(j)], [y(i), y(j)], 'k-', 'LineWidth', 1);
                end
            end
        end
    end
    
    % ================== 标注晶格矢量 ==================
    % 绘制a1基矢
    quiver(0, 0, a1(1), a1(2), 0, ...
           'k-', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(a1(1)/2, a1(2)/2, 'a_1', ...
         'FontSize', 12, 'FontWeight', 'bold');
    
    % 绘制a2基矢
    quiver(0, 0, a2(1), a2(2), 0, ...
           'k-', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(a2(1)/2, a2(2)/2, 'a_2', ...
         'FontSize', 12, 'FontWeight', 'bold');
    
    % ================== 设置图形属性 ==================
    title(sprintf('%dx%d 石墨烯超胞结构 (%d 原子)', Nx, Ny, size(atom_position, 1)));
    xlabel('x (Å)');
    ylabel('y (Å)');
    axis equal;  % 保持坐标轴比例一致
    grid on;     % 显示网格
    
    % 添加图例
    legend({'超胞边界', 'A子晶格', 'B子晶格'}, 'Location', 'best');
    
    % 设置坐标轴范围
    padding = a/2;  % 边界填充
    xlim([-padding, max(x) + padding]);
    ylim([-padding, max(y) + padding]);
    
    hold off;  % 结束绘图
end
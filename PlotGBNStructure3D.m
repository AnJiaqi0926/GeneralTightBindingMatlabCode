function PlotGBNStructure3D(a1, a2, a, atom_position, Nx, Ny)
    % 绘制石墨烯-氮化硼双层结构3D图
    %
    % 输入参数:
    %   a1, a2 - 超胞基矢 (2D)
    %   a - 晶格常数
    %   atom_position - 原子位置矩阵 [N×9]
    %   Nx, Ny - 超胞尺寸
    
    % 创建新图形窗口
    figure('Name', sprintf('%dx%d G/hBN双层结构', Nx, Ny), ...
           'Position', [100, 100, 1000, 800]);
    
    % 提取原子位置信息
    x = atom_position(:, 1);          % x坐标
    y = atom_position(:, 2);          % y坐标
    z = atom_position(:, 3);          % z坐标
    layer = atom_position(:, 4);      % 层数 (1=石墨烯, 2=hBN)
    sublattice = atom_position(:, 5); % 子晶格标识 
    atom_type = atom_position(:, 6);  % 原子类型 (1=碳, 3=氮, 4=硼)
    
    % 分离不同层和不同类型的原子
    % 石墨烯层 (layer=1)
    grapheneA_idx = (layer == 1) & (sublattice == 1);  % 上层石墨烯A子晶格
    grapheneB_idx = (layer == 1) & (sublattice == 2);  % 上层石墨烯B子晶格
    
    % hBN层 (layer=2)
    hBN_N_idx = (layer == 2) & (atom_type == 3);  % 下层石墨烯A子晶格
    hBN_B_idx = (layer == 2) & (atom_type == 4);  % 下层石墨烯B子晶格
    
    % 绘制图形
    hold on;
    
    % ================== 绘制超胞边界 ==================
    % 计算超胞顶点坐标 (底部和顶部)
    vertices_bottom = [0, 0, min(z); 
                      a1(1), a1(2), min(z);
                      a1(1) + a2(1), a1(2) + a2(2), min(z);
                      a2(1), a2(2), min(z);
                      0, 0, min(z)];
                  
    vertices_top = [0, 0, max(z); 
                  a1(1), a1(2), max(z);
                  a1(1) + a2(1), a1(2) + a2(2), max(z);
                  a2(1), a2(2), max(z);
                  0, 0, max(z)];
    
    % 绘制超胞边界 (底部和顶部)
    plot3(vertices_bottom(:, 1), vertices_bottom(:, 2), vertices_bottom(:, 3), ...
          'k-', 'LineWidth', 1.5);
    plot3(vertices_top(:, 1), vertices_top(:, 2), vertices_top(:, 3), ...
          'k-', 'LineWidth', 1.5);
    
    % 绘制超胞侧面连接线
    for i = 1:4
        plot3([vertices_bottom(i, 1), vertices_top(i, 1)], ...
              [vertices_bottom(i, 2), vertices_top(i, 2)], ...
              [vertices_bottom(i, 3), vertices_top(i, 3)], ...
              'k--', 'LineWidth', 1);
    end
    
    % ================== 绘制原子 ==================
    % 绘制石墨烯A原子 (红色)
    scatter3(x(grapheneA_idx), y(grapheneA_idx), z(grapheneA_idx), 100, ...
            'filled', 'MarkerFaceColor', [0.8, 0.2, 0.2], ...  % 红色
            'MarkerEdgeColor', 'k');
    
    % 绘制石墨烯B原子 (蓝色)
    scatter3(x(grapheneB_idx), y(grapheneB_idx), z(grapheneB_idx), 100, ...
            'filled', 'MarkerFaceColor', [0.2, 0.2, 0.8], ...  % 蓝色
            'MarkerEdgeColor', 'k');
    
    % 绘制hBN氮原子 (绿色)
    scatter3(x(hBN_N_idx), y(hBN_N_idx), z(hBN_N_idx), 100, ...
            'filled', 'MarkerFaceColor', [0.2, 0.8, 0.2], ...  % 绿色
            'MarkerEdgeColor', 'k');
    
    % 绘制hBN硼原子 (黄色)
    scatter3(x(hBN_B_idx), y(hBN_B_idx), z(hBN_B_idx), 100, ...
            'filled', 'MarkerFaceColor', [0.8, 0.8, 0.2], ...  % 黄色
            'MarkerEdgeColor', 'k');
    
    % ================== 绘制原子键 ==================
    % 仅对小型超胞绘制原子键 (避免性能问题)
    if Nx <= 10 && Ny <= 10
        d_nn = norm(a/sqrt(3));  % 最近邻距离
        
        % 遍历所有原子对
        for i = 1:size(atom_position, 1)
            for j = i+1:size(atom_position, 1)
                % 只绘制同一层内的键
                if layer(i) == layer(j)
                    % 计算原子间距离
                    dx = x(j) - x(i);
                    dy = y(j) - y(i);
                    dz = z(j) - z(i);
                    dist = sqrt(dx^2 + dy^2 + dz^2);
                    
                    % 如果是最近邻原子对，绘制键
                    if abs(dist - d_nn) < 0.1  % 允许0.1Å的容差
                        % 根据层数设置键的颜色
                        if layer(i) == 1  % 石墨烯层
                            line_color = [0.5, 0.5, 0.5];  % 灰色
                        else  % hBN层
                            line_color = [0.7, 0.7, 0.7];  % 浅灰色
                        end
                        
                        plot3([x(i), x(j)], [y(i), y(j)], [z(i), z(j)], ...
                              '-', 'Color', line_color, 'LineWidth', 1);
                    end
                end
            end
        end
    end
    
    % ================== 标注晶格矢量 ==================
    % 绘制a1基矢 (底部)
    quiver3(0, 0, min(z), a1(1), a1(2), 0, ...
            'k-', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    text(a1(1)/2, a1(2)/2, min(z), 'a_1', ...
         'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    
    % 绘制a2基矢 (底部)
    quiver3(0, 0, min(z), a2(1), a2(2), 0, ...
            'k-', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    text(a2(1)/2, a2(2)/2, min(z), 'a_2', ...
         'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    
    % 绘制z方向矢量
    quiver3(0, 0, min(z), 0, 0, max(z)-min(z), ...
            'r-', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    text(0, 0, (min(z)+max(z))/2, 'z', ...
         'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    
    % ================== 设置图形属性 ==================
    title(sprintf('%dx%d G/hBN双层结构 (%d 原子)', Nx, Ny, size(atom_position, 1)));
    xlabel('x (Å)');
    ylabel('y (Å)');
    zlabel('z (Å)');
    axis equal;  % 保持坐标轴比例一致
    grid on;     % 显示网格
    
    % 添加图例
    legend_items = {'超胞边界', '石墨烯-A', '石墨烯-B', '氮原子', '硼原子'};
    legend(legend_items, 'Location', 'best');
    
    % 设置坐标轴范围
    padding = a/2;  % 边界填充
    xlim([-padding, max(x) + padding]);
    ylim([-padding, max(y) + padding]);
    zlim([min(z)-1, max(z)+1]);  % z方向留出额外空间
    
    % 设置视角
    view(45, 30);  % 3D视角
    rotate3d on;   % 启用3D旋转
    
    hold off;  % 结束绘图
end
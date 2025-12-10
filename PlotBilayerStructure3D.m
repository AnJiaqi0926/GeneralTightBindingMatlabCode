function PlotBilayerStructure3D(a1, a2, a, atom_position, Nx, Ny, mode)
    % 绘制石墨烯-氮化硼双层结构3D图（支持多种显示模式）
    %
    % 输入参数:
    %   a1, a2 - 超胞基矢 (2D)
    %   a - 晶格常数
    %   atom_position - 原子位置矩阵 [N×9]
    %   Nx, Ny - 超胞尺寸
    %   mode - 显示模式: 
    %          'single' - 仅显示一个原胞
    %          '3cells' - 显示3个原胞
    %          '9cells' - 显示9个原胞
    
    % 创建新图形窗口
    figure('Name', sprintf('%dx%d G/hBN双层结构 (%s模式)', Nx, Ny, mode), ...
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
    grapheneA_idx = (sublattice == 1);  % 上层石墨烯A子晶格
    grapheneB_idx = (sublattice == 2);  % 上层石墨烯B子晶格
    
    % 绘制图形
    hold on;
    
    % ================== 根据模式确定要绘制的原胞 ==================
    if strcmp(mode, 'single')
        translations = [0, 0];
    elseif strcmp(mode, '3cells')
        translations = [-1, 0; 0, 0; 1, 0];
    elseif strcmp(mode, '9cells')
        translations = [-1, -1; -1, 0; -1, 1;
                         0, -1;  0, 0;  0, 1;
                         1, -1;  1, 0;  1, 1];
    else
        error('未知模式: %s。请使用 ''single'', ''3cells'' 或 ''9cells''。', mode);
    end
    
    % 为不同原胞设置不同的透明度
    if strcmp(mode, 'single')
        alphas = 1.0;
    elseif strcmp(mode, '3cells')
        alphas = [0.4, 1.0, 0.4];
    else % 9cells
        alphas = [0.1, 0.2, 0.1;
                  0.2, 1.0, 0.2;
                  0.1, 0.2, 0.1];
    end
    
    
    % ================== 绘制原子 ==================
    for idx = 1:size(translations, 1)
        i = translations(idx, 1);
        j = translations(idx, 2);
        
        % 计算平移矢量 (只在xy平面平移，z坐标不变)
        translation_vector = i * a1 + j * a2;
        
        % 计算平移后的原子位置
        x_translated = x + translation_vector(1);
        y_translated = y + translation_vector(2);
        z_translated = z; % z坐标不变
        
        % 获取当前原胞的透明度
        if strcmp(mode, 'single')
            current_alpha = alphas;
        elseif strcmp(mode, '3cells')
            current_alpha = alphas(idx);
        else % 9cells
            current_alpha = alphas(i+2, j+2); % +2是因为索引从1开始
        end
        
        % 绘制石墨烯A原子 (红色)
        scatter3(x_translated(grapheneA_idx), y_translated(grapheneA_idx), z_translated(grapheneA_idx), 100, ...
                'filled', 'MarkerFaceColor', [0.8, 0.2, 0.2], ...  % 红色
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', current_alpha, ...
                'MarkerEdgeAlpha', current_alpha*0.8);
        
        % 绘制石墨烯B原子 (蓝色)
        scatter3(x_translated(grapheneB_idx), y_translated(grapheneB_idx), z_translated(grapheneB_idx), 100, ...
                'filled', 'MarkerFaceColor', [0.2, 0.2, 0.8], ...  % 蓝色
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', current_alpha, ...
                'MarkerEdgeAlpha', current_alpha*0.8);
    end
    
    % ================== 绘制原子键 ==================
    % 仅对小型超胞绘制原子键 (避免性能问题)
    if Nx <= 10 && Ny <= 10
        d_nn = norm(a/sqrt(3));  % 最近邻距离
        
        % 收集所有要绘制的原子位置
        all_x = [];
        all_y = [];
        all_z = [];
        all_layer = [];
        
        for idx = 1:size(translations, 1)
            i = translations(idx, 1);
            j = translations(idx, 2);
            
            % 计算平移矢量
            translation_vector = i * a1 + j * a2;
            
            % 计算平移后的原子位置
            x_translated = x + translation_vector(1);
            y_translated = y + translation_vector(2);
            z_translated = z; % z坐标不变
            
            all_x = [all_x; x_translated];
            all_y = [all_y; y_translated];
            all_z = [all_z; z_translated];
            all_layer = [all_layer; layer];
        end
        
        % 绘制所有原子间的键
        for i = 1:length(all_x)
            for j = i+1:length(all_x)
                % 只绘制同一层内的键
                if all_layer(i) == all_layer(j)
                    % 计算原子间距离
                    dx = all_x(j) - all_x(i);
                    dy = all_y(j) - all_y(i);
                    dz = all_z(j) - all_z(i);
                    dist = sqrt(dx^2 + dy^2 + dz^2);
                    
                    % 如果是最近邻原子对，绘制键
                    if abs(dist - d_nn) < 0.1  % 允许0.1Å的容差
                        % 根据层数设置键的颜色
                        if all_layer(i) == 1  % 石墨烯层
                            line_color = [0.5, 0.5, 0.5];  % 灰色
                        else  % hBN层
                            line_color = [0.7, 0.7, 0.7];  % 浅灰色
                        end
                        
                        plot3([all_x(i), all_x(j)], [all_y(i), all_y(j)], [all_z(i), all_z(j)], ...
                              '-', 'Color', line_color, 'LineWidth', 1);
                    end
                end
            end
        end
    end
    
    % ================== 绘制超胞边界 ==================
    % 计算原胞顶点坐标 (底部和顶部)
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
    
    % ================== 绘制中心原胞边界 ==================
    % 绘制中心原胞边界（用粗线突出显示）
    plot3(vertices_bottom(:, 1), vertices_bottom(:, 2), vertices_bottom(:, 3), ...
          'r-', 'LineWidth', 3);
    plot3(vertices_top(:, 1), vertices_top(:, 2), vertices_top(:, 3), ...
          'r-', 'LineWidth', 3);
    
    % 绘制中心原胞侧面连接线
    for k = 1:4
        plot3([vertices_bottom(k, 1), vertices_top(k, 1)], ...
              [vertices_bottom(k, 2), vertices_top(k, 2)], ...
              [vertices_bottom(k, 3), vertices_top(k, 3)], ...
              'r-', 'LineWidth', 2);
    end
    
    % ================== 绘制周围原胞边界 ==================
    % 绘制周围原胞的边界（用虚线表示）
    for idx = 1:size(translations, 1)
        i = translations(idx, 1);
        j = translations(idx, 2);
        
        % 跳过中心原胞（已经在上面绘制了）
        if i == 0 && j == 0
            continue;
        end
        
        % 计算平移矢量
        translation_vector = i * a1 + j * a2;
        
        % 计算平移后的顶点坐标
        vertices_bottom_translated = vertices_bottom + [translation_vector, 0];
        vertices_top_translated = vertices_top + [translation_vector, 0];
        
        % 绘制边界（灰色虚线）
        plot3(vertices_bottom_translated(:, 1), vertices_bottom_translated(:, 2), vertices_bottom_translated(:, 3), ...
              'k--', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5, 0.6]);
        plot3(vertices_top_translated(:, 1), vertices_top_translated(:, 2), vertices_top_translated(:, 3), ...
              'k--', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5, 0.6]);
        
        % 绘制侧面连接线（虚线）
        for k = 1:4
            plot3([vertices_bottom_translated(k, 1), vertices_top_translated(k, 1)], ...
                  [vertices_bottom_translated(k, 2), vertices_top_translated(k, 2)], ...
                  [vertices_bottom_translated(k, 3), vertices_top_translated(k, 3)], ...
                  'k--', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5, 0.6]);
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
    title(sprintf('%dx%d G/hBN双层结构 (%s模式, 共%d原子)', Nx, Ny, mode, ...
          size(atom_position, 1)*size(translations, 1)));
    xlabel('x (Å)');
    ylabel('y (Å)');
    zlabel('z (Å)');
    axis equal;  % 保持坐标轴比例一致
    grid on;     % 显示网格
    
    % 添加图例
    legend_items = {'中心原胞边界', '石墨烯-A', '石墨烯-B', '氮原子', '硼原子', '邻胞边界'};
    legend(legend_items, 'Location', 'best');
    
    % 设置坐标轴范围
    padding = a/2;  % 边界填充
    x_center = mean(x);
    y_center = mean(y);
    x_range = max(x) - min(x);
    y_range = max(y) - min(y);
    
    % 根据模式调整坐标轴范围
    if strcmp(mode, 'single')
        xlim([-padding, max(x) + padding]);
        ylim([-padding, max(y) + padding]);
    elseif strcmp(mode, '3cells')
        xlim([x_center - 1.5*x_range - padding, x_center + 1.5*x_range + padding]);
        ylim([y_center - 0.5*y_range - padding, y_center + 0.5*y_range + padding]);
    else % 9cells
        xlim([x_center - 1.5*x_range - padding, x_center + 1.5*x_range + padding]);
        ylim([y_center - 1.5*y_range - padding, y_center + 1.5*y_range + padding]);
    end
    zlim([min(z)-1, max(z)+1]);  % z方向留出额外空间
    
    % 设置视角
    view(45, 30);  % 3D视角
    rotate3d on;   % 启用3D旋转
    
    hold off;  % 结束绘图
end
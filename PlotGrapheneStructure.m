function PlotGrapheneStructure(a1, a2, a, atom_position, Nx, Ny, mode)
    % 绘制石墨烯结构图函数（支持多种显示模式）
    %
    % 输入参数:
    %   a1, a2 - 超胞基矢
    %   a - 晶格常数
    %   atom_position - 原子位置矩阵
    %   Nx, Ny - 超胞尺寸
    %   mode - 显示模式: 
    %          'single' - 仅显示一个原胞
    %          '3cells' - 显示3个原胞
    %          '9cells' - 显示9个原胞
    % 创建新图形窗口
    figure('Name', sprintf('%dx%d 石墨烯超胞结构 (%s模式)', Nx, Ny, mode), ...
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
        alphas = [0.4, 0.6, 0.4;
                  0.6, 1.0, 0.6;
                  0.4, 0.6, 0.4];
    end
    
    % ================== 绘制原子 ==================
    for idx = 1:size(translations, 1)
        i = translations(idx, 1);
        j = translations(idx, 2);
        
        % 计算平移矢量
        translation_vector = i * a1 + j * a2;
        
        % 计算平移后的原子位置
        x_translated = x + translation_vector(1);
        y_translated = y + translation_vector(2);
        
        % 获取当前原胞的透明度
        if strcmp(mode, 'single')
            current_alpha = alphas;
        elseif strcmp(mode, '3cells')
            current_alpha = alphas(idx);
        else % 9cells
            current_alpha = alphas(i+2, j+2); % +2是因为索引从1开始
        end
        
        % 绘制A子晶格原子
        scatter(x_translated(idxA), y_translated(idxA), 80, 'filled', ...
                'MarkerFaceColor', [0.8, 0.2, 0.2], ... % 红色
                'MarkerEdgeColor', [0.5, 0.1, 0.1], ...
                'MarkerFaceAlpha', current_alpha, ...
                'MarkerEdgeAlpha', current_alpha*0.8);
        
        % 绘制B子晶格原子
        scatter(x_translated(idxB), y_translated(idxB), 80, 'filled', ...
                'MarkerFaceColor', [0.2, 0.2, 0.8], ... % 蓝色
                'MarkerEdgeColor', [0.1, 0.1, 0.5], ...
                'MarkerFaceAlpha', current_alpha, ...
                'MarkerEdgeAlpha', current_alpha*0.8);
    end
    
    % ================== 绘制原子键（全局） ==================
    % 仅对小型超胞绘制原子键 (避免性能问题)
    if Nx <= 10 && Ny <= 10
        d_nn = norm(a/sqrt(3));  % 最近邻距离
        
        % 收集所有要绘制的原子位置
        all_x = [];
        all_y = [];
        
        for idx = 1:size(translations, 1)
            i = translations(idx, 1);
            j = translations(idx, 2);
            
            % 计算平移矢量
            translation_vector = i * a1 + j * a2;
            
            % 计算平移后的原子位置
            x_translated = x + translation_vector(1);
            y_translated = y + translation_vector(2);
            
            all_x = [all_x; x_translated];
            all_y = [all_y; y_translated];
        end
        
        % 绘制所有原子间的键
        for i = 1:length(all_x)
            for j = i+1:length(all_x)
                % 计算原子间距离
                dx = all_x(j) - all_x(i);
                dy = all_y(j) - all_y(i);
                dist = sqrt(dx^2 + dy^2);
                
                % 如果是最近邻原子对，绘制键
                if abs(dist - d_nn) < 0.1  % 允许0.1Å的容差
                    plot([all_x(i), all_x(j)], [all_y(i), all_y(j)], 'k-', 'LineWidth', 1.5);
                end
            end
        end
    end
    
    % ================== 绘制原胞边界 ==================
    % 计算原胞顶点坐标
    vertices = [0, 0; 
                a1(1), a1(2);
                a1(1) + a2(1), a1(2) + a2(2);
                a2(1), a2(2);
                0, 0];
    
    % 绘制中心原胞边界（用粗红线突出显示）
    plot(vertices(:, 1), vertices(:, 2), 'r-', 'LineWidth', 3);
    
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
        vertices_translated = vertices + translation_vector;
        
        % 绘制边界（灰色虚线）
        plot(vertices_translated(:, 1), vertices_translated(:, 2), 'k--', ...
             'LineWidth', 1, 'Color', [0.5, 0.5, 0.5, 0.6]);
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
    title(sprintf('%dx%d 石墨烯超胞结构 (%s模式, 共%d原子)', Nx, Ny, mode, ...
          size(atom_position, 1)*size(translations, 1)));
    xlabel('x (Å)');
    ylabel('y (Å)');
    axis equal;  % 保持坐标轴比例一致
    grid on;     % 显示网格
    
    % 添加图例
    legend({'中心原胞边界', 'A子晶格', 'B子晶格', '邻胞边界'}, ...
           'Location', 'best');
    
    % 设置坐标轴范围（扩展到显示所有原胞）
    padding = a;  % 增加边界填充
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
    
    hold off;  % 结束绘图
end
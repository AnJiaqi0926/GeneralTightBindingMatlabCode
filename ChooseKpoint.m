function [b1, b2, kx_v, ky_v, kpath, high_sym_points] = ChooseKpoint(a1, a2, nsum, mode, custom_high_sym_points)
    % 选择k点路径函数
    %
    % 输入参数:
    %   a1, a2 - 正格子基矢
    %   nsum   - k点数量
    %   mode   - 'full' (全布里渊区) 或 'high_sym' (高对称线)
    %   custom_high_sym_points - (可选) 自定义高对称点路径 [n×2矩阵]
    %
    % 输出参数:
    %   b1, b2 - 倒格子基矢
    %   kx_v, ky_v - k点笛卡尔坐标
    %   kpath - 路径长度坐标 (仅用于高对称线模式)
    %   high_sym_points - 使用的高对称点分数坐标
    
    % 设置默认模式
    if nargin < 4
        mode = 'full';  % 默认使用全布里渊区模式
    end
    
    % ================== 计算倒格矢 ==================
    % 构建基矢矩阵
    A = [a1(:), a2(:)];
    
    % 检查基矢是否线性无关
    det_A = det(A);
    if abs(det_A) < eps
        error('基矢a1和a2线性相关，无法构成二维晶格。');
    end
    
    % 计算倒格矢矩阵: B = 2π * inv(A)^T
    inv_A_T = inv(A)';
    B = 2 * pi * inv_A_T;
    
    % 提取倒格矢
    b1 = B(:, 1).';  % 第一倒格矢
    b2 = B(:, 2).';  % 第二倒格矢
    
    % ================== 高对称点定义 ==================
    % 默认高对称点路径 (六角晶格)
    default_high_sym_points = [0, 0;       % Γ点
                               0, 1/2;     % M点  
                               1/3, 2/3;   % K点
                               0, 0];      % Γ点 (闭合路径)
    
    % 根据输入选择高对称点路径
    if nargin >= 5 && ~isempty(custom_high_sym_points)
        high_sym_points = custom_high_sym_points;
        fprintf('使用自定义高对称点路径...\n');
    else
        high_sym_points = default_high_sym_points;
        fprintf('使用默认高对称点路径 (Γ-M-K-Γ)...\n');
    end
    
    % 检查高对称点路径的有效性
    if strcmp(mode, 'high_sym') && size(high_sym_points, 1) < 2
        error('高对称点路径至少需要2个点');
    end
    
    % ================== k点路径生成 ==================
    if strcmp(mode, 'high_sym')
        % 高对称线模式
        n_segments = size(high_sym_points, 1) - 1;  % 路径段数
        
        % 计算每段k点数 (根据路径长度分配)
        segment_lengths = zeros(1, n_segments);
        total_length = 0;
        
        % 预计算各段长度
        for seg_idx = 1:n_segments
            start_point = high_sym_points(seg_idx, :);
            end_point = high_sym_points(seg_idx + 1, :);
            
            % 转换为笛卡尔坐标计算距离
            start_cart = start_point(1) * b1 + start_point(2) * b2;
            end_cart = end_point(1) * b1 + end_point(2) * b2;
            
            segment_lengths(seg_idx) = norm(end_cart - start_cart);
            total_length = total_length + segment_lengths(seg_idx);
        end
        
        % 根据长度比例分配k点数
        points_per_segment = zeros(1, n_segments);
        for seg_idx = 1:n_segments
            points_per_segment(seg_idx) = max(2, round(nsum * segment_lengths(seg_idx) / total_length));
        end
        
        % 调整总点数
        total_points = sum(points_per_segment);
        if total_points > nsum
            % 按比例减少
            scale_factor = nsum / total_points;
            points_per_segment = max(2, round(points_per_segment * scale_factor));
        end
        
        kx_v = [];
        ky_v = [];
        kpath = [];
        cumulative_length = 0;
        
        % 生成各路径段的k点
        for seg_idx = 1:n_segments
            start_point = high_sym_points(seg_idx, :);
            end_point = high_sym_points(seg_idx + 1, :);
            
            % 在分数坐标中线性插值
            frac_kx = linspace(start_point(1), end_point(1), points_per_segment(seg_idx));
            frac_ky = linspace(start_point(2), end_point(2), points_per_segment(seg_idx));
            
            % 转换为笛卡尔坐标
            segment_kx = frac_kx * b1(1) + frac_ky * b2(1);
            segment_ky = frac_kx * b1(2) + frac_ky * b2(2);
            
            % 计算当前路径段的长度
            segment_length = norm([segment_kx(end) - segment_kx(1), segment_ky(end) - segment_ky(1)]);
            
            % 生成路径长度坐标
            kpath_segment = linspace(cumulative_length, cumulative_length + segment_length, points_per_segment(seg_idx));
            
            % 添加到总路径 (避免重复点)
            if seg_idx > 1
                kx_v = [kx_v, segment_kx(2:end)];
                ky_v = [ky_v, segment_ky(2:end)];
                kpath = [kpath, kpath_segment(2:end)];
            else
                kx_v = [kx_v, segment_kx];
                ky_v = [ky_v, segment_ky];
                kpath = [kpath, kpath_segment];
            end
            
            cumulative_length = cumulative_length + segment_length;
        end
        
    else
        fprintf('生成全布里渊区k点网格...\n');
        
        kpath = [];
        
        % 计算布里渊区边界
        % 使用倒格子基矢定义布里渊区的边界
        % 对于任意晶格，我们生成覆盖整个倒格子原胞的k点
        
        % 确定网格尺寸
        nside = nsum;
        
        % 在分数坐标中生成均匀网格
        % 范围从-0.5到0.5，覆盖整个倒格子原胞
        u = linspace(-0.5, 0.5, nside);
        v = linspace(-0.5, 0.5, nside);
        
        [U, V] = meshgrid(u, v);
        u_vec = U(:);
        v_vec = V(:);
        
        % 转换为笛卡尔坐标: k = u*b1 + v*b2
        kx_v = u_vec * b1(1) + v_vec * b2(1);
        ky_v = u_vec * b1(2) + v_vec * b2(2);
        
        kx_v = kx_v';
        ky_v = ky_v';
    end
end
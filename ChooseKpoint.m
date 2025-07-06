function [b1, b2, kx_v, ky_v, kpath, high_sym_points] = ChooseKpoint(a1, a2, nsum, mode)
    % 选择k点路径函数
    %
    % 输入参数:
    %   a1, a2 - 正格子基矢
    %   nsum   - k点数量
    %   mode   - 'full' (全布里渊区) 或 'high_sym' (高对称线)
    %
    % 输出参数:
    %   b1, b2 - 倒格子基矢
    %   kx_v, ky_v - k点笛卡尔坐标
    %   kpath - 路径长度坐标 (仅用于高对称线模式)
    %   high_sym_points - 高对称点分数坐标
    
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
    % 分数坐标 (以正格子基矢为基准)
    Gamma = [0, 0];    % Γ点
    M = [0, 1/2];    % M点
    K = [1/3, 2/3];    % K点
    high_sym_points = [Gamma; M; K; Gamma];  % Γ -> M -> K -> Γ 路径
    
    % ================== k点路径生成 ==================
    if strcmp(mode, 'high_sym')
        % 高对称线模式: 沿Γ-M-K-Γ路径生成k点
        fprintf('生成高对称线k点路径 (Γ-M-K-Γ)...\n');
        
        n_segments = size(high_sym_points, 1) - 1;  % 路径段数
        points_per_segment = ceil(nsum / n_segments);  % 每段k点数
        
        kx_v = [];  % 初始化k_x坐标
        ky_v = [];  % 初始化k_y坐标
        kpath = []; % 初始化路径长度
        cumulative_length = 0;  % 累计路径长度
        
        % 遍历每条路径段
        for seg_idx = 1:n_segments
            start_point = high_sym_points(seg_idx, :);      % 起点 (分数坐标)
            end_point = high_sym_points(seg_idx + 1, :);    % 终点 (分数坐标)
            
            % 在分数坐标中线性插值
            frac_kx = linspace(start_point(1), end_point(1), points_per_segment);
            frac_ky = linspace(start_point(2), end_point(2), points_per_segment);
            
            % 转换为笛卡尔坐标: k_cart = kx_frac * b1 + ky_frac * b2
            segment_kx = frac_kx * b1(1) + frac_ky * b2(1);
            segment_ky = frac_kx * b1(2) + frac_ky * b2(2);
            
            % 计算当前路径段的长度
            segment_length = sqrt((segment_kx(end) - segment_kx(1))^2 + ...
                                  (segment_ky(end) - segment_ky(1))^2);
            
            % 生成路径长度坐标
            kpath_segment = linspace(cumulative_length, ...
                                    cumulative_length + segment_length, ...
                                    points_per_segment);
            
            % 添加到总路径
            kx_v = [kx_v, segment_kx];
            ky_v = [ky_v, segment_ky];
            kpath = [kpath, kpath_segment];
            
            % 更新累计路径长度
            cumulative_length = cumulative_length + segment_length;
        end
        
        % 确保总点数不超过nsum
        if length(kx_v) > nsum
            kx_v = kx_v(1:nsum);
            ky_v = ky_v(1:nsum);
            kpath = kpath(1:nsum);
        end
        
    else
        % 全布里渊区模式: 生成均匀网格k点
        fprintf('生成全布里渊区k点网格...\n');
        
        kpath = [];  % 路径长度坐标不适用
        
        % 在倒格子空间生成均匀网格
        kx_v = linspace(-b1(1), b1(1), nsum);
        ky_v = linspace(-b2(2), b2(2), nsum);
        
        % 创建网格并向量化
        [KX, KY] = meshgrid(kx_v, ky_v);
        kx_v = KX(:)';
        ky_v = KY(:)';
    end
end
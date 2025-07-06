function [b1, b2, k_v] = ChooseKpoint_Ribbon(a1, a2, nsum)
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
    
    % ================== k点路径生成 ==================
    % 高对称线模式: 沿Γ-M-K-Γ路径生成k点
    fprintf('生成高对称线k点路径 (-K --- K)...\n');    
    k_v = linspace(0,2*pi/b1(1),nsum);  % 初始化k_x坐标
end
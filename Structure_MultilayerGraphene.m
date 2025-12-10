function [a1, a2, atom_position] = Structure_MultilayerGraphene(Nx, Ny, plot_structure,params)
    % 构建石墨烯超胞函数
    %
    % 输入参数:
    %   Nx, Ny - 超胞在a1和a2方向上的重复次数
    %   plot_structure - 是否绘制原子结构图
    %
    % 输出参数:
    %   a1, a2 - 超胞基矢
    %   atom_position - 原子位置矩阵 [N×9]
    
    % ================== 参数设置 ==================
    if isfield(params, 'ifspin')
        ifspin = params.ifspin;
    else
        ifspin = false;
    end

    if isfield(params, 'bulkcellMode')
        bulkcellMode = params.bulkcellMode;
    else
        bulkcellMode = "9cells";
    end

    if isfield(params, 'a')
        a = params.a;
    else
        a = 2.46;
    end

    if isfield(params, 'StackedMode')
        StackedMode = params.StackedMode;
    else
        StackedMode = "AB";
    end

    if isfield(params, 'total_layer')
        total_layer = params.total_layer;
    else
        total_layer = 4;
    end
    
    % ================== 函数本体 ==================    
    % 石墨烯原胞基矢 (正交表示)
    a1_single = [a, 0];                % a1基矢
    a2_single = [a/2, a*sqrt(3)/2];    % a2基矢
    
    % 超胞基矢
    a1 = Nx * a1_single;  % 沿a1方向放大Nx倍
    a2 = Ny * a2_single;  % 沿a2方向放大Ny倍
    
    % 原胞内原子位置 (分数坐标)
    % A原子位置: (1/3, 1/3)
    % B原子位置: (2/3, 2/3)
    
    posA = (a1_single + a2_single)/3;
    posB = 2*(a1_single + a2_single)/3; 
    shift = (a1_single + a2_single)/3;

    % ================== 构建超胞 ==================
    if ifspin
        total_atoms = 4 * total_layer * Nx * Ny;  % 总原子数 (每个原胞4个原子)
        fprintf('构建 %dx%d 石墨烯超胞: %d 个原子\n', Nx, Ny, total_atoms);
    else
        total_atoms = 2 * total_layer * Nx * Ny;  % 总原子数 (每个原胞4个原子)
        fprintf('构建 %dx%d 石墨烯超胞: %d 个原子\n', Nx, Ny, total_atoms);
    end
        
    % 预分配原子位置矩阵
    atom_position = zeros(total_atoms, 10);
    atom_count = 0;  % 原子计数器
    
    % 在二维网格上放置原子
    for ix = 0:Nx-1       % 沿a1方向
        for iy = 0:Ny-1   % 沿a2方向
            % 计算当前原胞的位移
            displacement = ix * a1_single + iy * a2_single;
            for iz = 0:total_layer-1
                layer_shift = getLayerShift(iz, StackedMode, shift);
                % 添加A原子 (子晶格A)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posA(1) + layer_shift(1), ... % x坐标
                    displacement(2) + posA(2) + layer_shift(2), ... % y坐标
                    3.35*iz, ...                   % z坐标
                    iz + 1, ...                    % 层数
                    1, ...                         % 子晶格 (A=1)
                    1, ...                         % 原子类型 (碳=1)
                    1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正 (未使用)
                
                % 添加B原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posB(1) + layer_shift(1), ... % x坐标
                    displacement(2) + posB(2) + layer_shift(2), ... % y坐标
                    3.35*iz, ...                   % z坐标
                    iz + 1, ...                    % 层数
                    2, ...                         % 子晶格 (B=2)
                    1, ...                         % 原子类型
                    1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正
    
                if ifspin
                    % 添加A原子 (子晶格A)
                    atom_count = atom_count + 1;
                    atom_position(atom_count, :) = [...
                        displacement(1) + posA(1) + layer_shift(1), ... % x坐标
                        displacement(2) + posA(2) + layer_shift(2), ... % y坐标
                        3.35*iz, ...                   % z坐标
                        iz + 1, ...                    % 层数
                        1, ...                         % 子晶格 (A=1)
                        1, ...                         % 原子类型 (碳=1)
                        -1, ...                        % 自旋类型 (1，-1)
                        0, 0, 0];                      % 位移修正 (未使用)
                    
                    % 添加B原子 (子晶格B)
                    atom_count = atom_count + 1;
                    atom_position(atom_count, :) = [...
                        displacement(1) + posB(1) + layer_shift(1), ... % x坐标
                        displacement(2) + posB(2) + layer_shift(2), ... % y坐标
                        3.35*iz, ...                   % z坐标
                        iz + 1, ...                    % 层数
                        2, ...                         % 子晶格 (B=2)
                        1, ...                         % 原子类型
                        -1, ...                        % 自旋类型 (1，-1)
                        0, 0, 0];                      % 位移修正
                end
            end
        end
    end
    
    % ================== 可视化结构 ==================
    if plot_structure
        PlotBilayerStructure3D(a1, a2, a, atom_position, Nx, Ny,bulkcellMode);
    end
end

% ================== 辅助函数：计算层偏移 ==================
function layer_shift = getLayerShift(iz, StackedMode, shift)
    % 根据堆叠模式和层索引计算层偏移
    %
    % 输入参数:
    %   iz - 层索引（从0开始）
    %   StackedMode - 堆叠模式 ('AA', 'AB', 'ABC')
    %   shift - 基本偏移向量 (a1+a2)/3
    %
    % 输出参数:
    %   layer_shift - 当前层的偏移向量
    
    switch StackedMode
        case 'AA'
            % AA堆叠：所有层原子垂直对齐，无偏移
            layer_shift = [0, 0];
            
        case 'AB'
            % AB堆叠：交替层偏移
            % 偶数层无偏移，奇数层有偏移
            mod2 = mod(iz, 2);
            layer_shift = mod2 * shift;       % 奇数层
            
        case 'ABC'
            % ABC堆叠：三层一个周期
            % 第0层: 无偏移
            % 第1层: shift
            % 第2层: 2*shift
            % 第3层: 无偏移 (因为3*shift = a1+a2, 等价于晶格矢量)
            mod3 = mod(iz, 3);
            layer_shift = mod3 * shift;
            
        otherwise
            error('不支持的堆叠模式: %s', StackedMode);
    end
end
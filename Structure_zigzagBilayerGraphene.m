function [a1, a2, atom_position] = Structure_zigzagBilayerGraphene(Nx, Ny, plot_structure,params)
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
    if nargin < 3
        plot_structure = true;  % 默认绘制结构图
    end

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
    
    % ================== 函数本体 ==================      
    % 石墨烯原胞基矢 (正交表示)
    a1_single = [a, 0];                % a1基矢
    a2_single = [0, a*sqrt(3)];    % a2基矢
    
    % 超胞基矢
    a1 = Nx * a1_single;  % 沿a1方向放大Nx倍
    a2 = Ny * a2_single;  % 沿a2方向放大Ny倍
    
    % 原胞内原子位置 (分数坐标)
    % A原子位置: (1/3, 1/3)
    % B原子位置: (2/3, 2/3)
    posA = [a/2,a/(2*sqrt(3))];
    posB = [a,a/sqrt(3)];
    posC = [a,2*a/sqrt(3)];
    posD = [a/2,5*a/(2*sqrt(3))];
    if strcmp(StackedMode, "AA")
        posE = [a/2,a/(2*sqrt(3))];
        posF = [a,a/sqrt(3)];
        posG = [a,2*a/sqrt(3)];
        posH = [a/2,5*a/(2*sqrt(3))];
    elseif strcmp(StackedMode, "AB")
        posE = [a/2,a/(2*sqrt(3))+a/sqrt(3)];
        posF = [a,a/sqrt(3)+a/sqrt(3)];
        posG = [a,2*a/sqrt(3)+a/sqrt(3)];
        posH = [a/2,5*a/(2*sqrt(3))+a/sqrt(3)];
    elseif strcmp(StackedMode, "BA")
        posE = [a/2,a/(2*sqrt(3))+2*a/sqrt(3)];
        posF = [a,a/sqrt(3)+2*a/sqrt(3)];
        posG = [a,2*a/sqrt(3)-a/sqrt(3)];
        posH = [a/2,5*a/(2*sqrt(3))-a/sqrt(3)];
    end
    
    % ================== 构建超胞 ==================
    if ifspin
        total_atoms = 8 * Nx * Ny;  % 总原子数 (每个原胞4个原子)
        fprintf('构建 %dx%d 石墨烯超胞: %d 个原子\n', Nx, Ny, total_atoms);
    else
        total_atoms = 4 * Nx * Ny;  % 总原子数 (每个原胞4个原子)
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
            
            % 添加A原子 (子晶格A)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posA(1), ... % x坐标
                displacement(2) + posA(2), ... % y坐标
                0, ...                         % z坐标 (二维材料设为0)
                1, ...                         % 层数 (单层设为1)
                1, ...                         % 子晶格 (A=1)
                1, ...                         % 原子类型 (碳=1)
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正 (未使用)
            
            % 添加B原子 (子晶格B)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posB(1), ... % x坐标
                displacement(2) + posB(2), ... % y坐标
                0, ...                         % z坐标
                1, ...                         % 层数
                2, ...                         % 子晶格 (B=2)
                1, ...                         % 原子类型
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正

            % 添加C原子 (子晶格A)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posC(1), ... % x坐标
                displacement(2) + posC(2), ... % y坐标
                0, ...                         % z坐标 (二维材料设为0)
                1, ...                         % 层数 (单层设为1)
                1, ...                         % 子晶格 (A=1)
                1, ...                         % 原子类型 (碳=1)
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正 (未使用)
            
            % 添加D原子 (子晶格B)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posD(1), ... % x坐标
                displacement(2) + posD(2), ... % y坐标
                0, ...                         % z坐标
                1, ...                         % 层数
                2, ...                         % 子晶格 (B=2)
                1, ...                         % 原子类型
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正

            % 添加A原子 (子晶格A)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posE(1), ... % x坐标
                displacement(2) + posE(2), ... % y坐标
                0, ...                         % z坐标 (二维材料设为0)
                1, ...                         % 层数 (单层设为1)
                1, ...                         % 子晶格 (A=1)
                1, ...                         % 原子类型 (碳=1)
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正 (未使用)
            
            % 添加B原子 (子晶格B)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posF(1), ... % x坐标
                displacement(2) + posF(2), ... % y坐标
                0, ...                         % z坐标
                1, ...                         % 层数
                2, ...                         % 子晶格 (B=2)
                1, ...                         % 原子类型
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正

            % 添加C原子 (子晶格A)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posG(1), ... % x坐标
                displacement(2) + posG(2), ... % y坐标
                0, ...                         % z坐标 (二维材料设为0)
                1, ...                         % 层数 (单层设为1)
                1, ...                         % 子晶格 (A=1)
                1, ...                         % 原子类型 (碳=1)
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正 (未使用)
            
            % 添加D原子 (子晶格B)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posH(1), ... % x坐标
                displacement(2) + posD(2), ... % y坐标
                0, ...                         % z坐标
                1, ...                         % 层数
                2, ...                         % 子晶格 (B=2)
                1, ...                         % 原子类型
                1, ...                         % 自旋类型 (1，-1)
                0, 0, 0];                      % 位移修正

            if ifspin
                % 添加A原子 (子晶格A)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posA(1), ... % x坐标
                    displacement(2) + posA(2), ... % y坐标
                    0, ...                         % z坐标 (二维材料设为0)
                    1, ...                         % 层数 (单层设为1)
                    1, ...                         % 子晶格 (A=1)
                    1, ...                         % 原子类型 (碳=1)
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正 (未使用)
                
                % 添加B原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posB(1), ... % x坐标
                    displacement(2) + posB(2), ... % y坐标
                    0, ...                         % z坐标
                    1, ...                         % 层数
                    2, ...                         % 子晶格 (B=2)
                    1, ...                         % 原子类型
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正
    
                % 添加C原子 (子晶格A)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posC(1), ... % x坐标
                    displacement(2) + posC(2), ... % y坐标
                    0, ...                         % z坐标 (二维材料设为0)
                    1, ...                         % 层数 (单层设为1)
                    1, ...                         % 子晶格 (A=1)
                    1, ...                         % 原子类型 (碳=1)
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正 (未使用)
                
                % 添加D原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posD(1), ... % x坐标
                    displacement(2) + posD(2), ... % y坐标
                    0, ...                         % z坐标
                    1, ...                         % 层数
                    2, ...                         % 子晶格 (B=2)
                    1, ...                         % 原子类型
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正

                % 添加A原子 (子晶格A)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posE(1), ... % x坐标
                    displacement(2) + posE(2), ... % y坐标
                    0, ...                         % z坐标 (二维材料设为0)
                    1, ...                         % 层数 (单层设为1)
                    1, ...                         % 子晶格 (A=1)
                    1, ...                         % 原子类型 (碳=1)
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正 (未使用)
                
                % 添加B原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posF(1), ... % x坐标
                    displacement(2) + posF(2), ... % y坐标
                    0, ...                         % z坐标
                    1, ...                         % 层数
                    2, ...                         % 子晶格 (B=2)
                    1, ...                         % 原子类型
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正
    
                % 添加C原子 (子晶格A)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posG(1), ... % x坐标
                    displacement(2) + posG(2), ... % y坐标
                    0, ...                         % z坐标 (二维材料设为0)
                    1, ...                         % 层数 (单层设为1)
                    1, ...                         % 子晶格 (A=1)
                    1, ...                         % 原子类型 (碳=1)
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正 (未使用)
                
                % 添加D原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posH(1), ... % x坐标
                    displacement(2) + posH(2), ... % y坐标
                    0, ...                         % z坐标
                    1, ...                         % 层数
                    2, ...                         % 子晶格 (B=2)
                    1, ...                         % 原子类型
                    -1, ...                         % 自旋类型 (1，-1)
                    0, 0, 0];                      % 位移修正
            end
        end
    end
    
    % ================== 可视化结构 ==================
    if plot_structure
        PlotGrapheneStructure(a1, a2, a, atom_position, Nx, Ny,"3cells");
    end
end
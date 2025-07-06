function [a1, a2, atom_position] = Structure_GBNMoirecell(Nx, Ny, plot_structure)
    % 构建石墨烯-氮化硼超胞函数
    %
    % 输入参数:
    %   Nx, Ny - 超胞在a1和a2方向上的重复次数
    %   plot_structure - 是否绘制原子结构图
    %
    % 输出参数:
    %   a1, a2 - 超胞基矢
    %   atom_position - 原子位置矩阵 [N×10]
    
    % ================== 参数设置 ==================
    if nargin < 3
        plot_structure = true;  % 默认绘制结构图
    end
    
    a = 2.46;  % 石墨烯晶格常数 (Å)
    aBN = a * Nx/(Nx-1); % 氮化硼晶格常数 (Å)
    
    % 石墨烯原胞基矢 (正交表示)
    a1_single = [a, 0];                % a1基矢
    a2_single = [a/2, a*sqrt(3)/2];    % a2基矢
    a1_BNsingle = [aBN, 0];                % a1基矢
    a2_BNsingle = [aBN/2, aBN*sqrt(3)/2];    % a2基矢
    
    % 超胞基矢
    a1 = Nx * a1_single;  % 沿a1方向放大Nx倍
    a2 = Ny * a2_single;  % 沿a2方向放大Ny倍
    
    % 原胞内原子位置 (分数坐标)
    % A原子位置: (1/3, 1/3)
    % B原子位置: (2/3, 2/3)
    posA = (a1_single + a2_single)/3;
    posB = 2*(a1_single + a2_single)/3;
    posC = (a1_BNsingle + a2_BNsingle)/3;
    posD = 2*(a1_BNsingle + a2_BNsingle)/3;
    
    % ================== 构建超胞 ==================
    total_atoms = 2 * Nx * Ny + 2 * (Nx-1) * (Ny-1);  % 总原子数 (每个原胞4个原子)
    fprintf('构建 %dx%d 石墨烯-氮化硼超胞: %d 个原子\n', Nx, Ny, total_atoms);
    
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
                50.0, ...                         % z坐标 (第一层50.0)
                1, ...                         % 层数 (第一层)
                1, ...                         % 子晶格 (A=1)
                1, ...                         % 原子类型 (碳=1)
                1, ...                         % 自旋类型 (1，-1)
                (displacement(1) + posA(1))/a1(1), (displacement(2) + posA(2))/a1(1), 0];                      % 位移修正 (未使用)
            
            % 添加B原子 (子晶格B)
            atom_count = atom_count + 1;
            atom_position(atom_count, :) = [...
                displacement(1) + posB(1), ... % x坐标
                displacement(2) + posB(2), ... % y坐标
                50.0, ...                         % z坐标 (第一层50.0) 
                1, ...                         % 层数 (第一层) 
                2, ...                         % 子晶格 (B=2)
                1, ...                         % 原子类型
                1, ...                         % 自旋类型 (1，-1)
                (displacement(1) + posB(1))/a1(1), (displacement(2) + posB(2))/a1(1), 0];                      % 位移修正

            if ix ~= Nx-1 && iy ~= Ny-1 
                displacement = ix * a1_BNsingle + iy * a2_BNsingle;
                % 添加C原子 (子晶格N)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posC(1), ... % x坐标
                    displacement(2) + posC(2), ... % y坐标
                    53.05, ...                         % z坐标 (第二层53.35)
                    2, ...                         % 层数(第二层)
                    3, ...                         % 子晶格 (N=3)
                    3, ...                         % 原子类型
                    1, ...                         % 自旋类型 (1，-1)
                    (displacement(1) + posC(1))/a1(1), (displacement(2) + posC(2))/a1(1), 0];                      % 位移修正
    
                % 添加D原子 (子晶格B)
                atom_count = atom_count + 1;
                atom_position(atom_count, :) = [...
                    displacement(1) + posD(1), ... % x坐标
                    displacement(2) + posD(2), ... % y坐标
                    53.05, ...                         % z坐标 (第二层53.35)
                    2, ...                         % 层数(第二层)
                    4, ...                         % 子晶格 (B=4)
                    4, ...                         % 原子类型
                    1, ...                         % 自旋类型 (1，-1)
                    (displacement(1) + posD(1))/a1(1), (displacement(2) + posD(2))/a1(1), 0];                      % 位移修正
            end
        end
    end
    
    % ================== 可视化结构 ==================
    if plot_structure
        PlotGBNStructure3D(a1, a2, a, atom_position, Nx, Ny);
    end
end
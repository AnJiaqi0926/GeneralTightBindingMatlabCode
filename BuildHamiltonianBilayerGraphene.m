function [Hamiltonian,neighbor_cell, disp_vec] = BuildHamiltonianBilayerGraphene(a1, a2, atom_position, N, neighbor_cell, params)
    % 构建紧束缚模型哈密顿量函数
    %
    % 输入参数:
    %   a1, a2 - 正格子基矢
    %   atom_position - 原子位置矩阵 [N×9]
    %   N - 原子总数
    %   t - 跃迁积分强度
    %
    % 输出参数:
    %   Hamiltonian - 哈密顿量张量 [N×N×N_neighbor]
    %   neighbor_cell - 近邻原胞列表
    %   disp_vec - 位移矢量
    
    % ================== 参数设置 ==================
    if isfield(params, 'a')
        a_lattice = params.a;
    else
        a_lattice = 2.46;
    end
    a_0 = a_lattice/sqrt(3);
    d_0 = 3.35;
    r0 = 0.184*a_lattice;
    

    % ================== S-K 参数 (GBN) ==================
    if isfield(params, 't1')
        t1 = params.t1;
    else
        t1 = -2.7;
    end

    if isfield(params, 't2')
        t2 = params.t2;
    else
        t2 = 0.48;
    end

    if isfield(params, 'E_elec')
        E_elec = params.E_elec;
    else
        E_elec = 0.0;
    end

    
    % ================== 近邻原胞定义 ==================
    % 定义9个最近邻原胞 (包括原胞自身)
    N_neighbor = size(neighbor_cell, 1);  % 近邻原胞数量
    
    % ================== 预计算位移矢量 ==================
    % 每个近邻原胞对应的位移矢量: R = n1*a1 + n2*a2
    disp_vec = zeros(N_neighbor, 2);  % 初始化位移矢量矩阵
    for k = 1:N_neighbor
        n1 = neighbor_cell(k, 1);  % a1方向倍数
        n2 = neighbor_cell(k, 2);  % a2方向倍数
        disp_vec(k, 1) = n1 * a1(1) + n2 * a2(1);  % x分量
        disp_vec(k, 2) = n1 * a1(2) + n2 * a2(2);  % y分量
    end
    
    % ================== 提取原子坐标 ==================
    % 预提取原子坐标以提高效率
    x = atom_position(:, 1);  % x坐标
    y = atom_position(:, 2);  % y坐标
    z = atom_position(:, 3);  % z坐标
    layer = atom_position(:, 4);  % 层数
    sublattice = atom_position(:, 5);  % 子晶格
    atom_type = atom_position(:, 6);  % 原子种类
    dx = atom_position(:, 7);  % 原子种类
    dy = atom_position(:, 8);  % 原子种类
    d = atom_position(:, 9);  % 原子种类

    
    % ================== 构建哈密顿量 ==================
    % 初始化哈密顿量张量
    Hamiltonian = zeros(N, N, N_neighbor);
    
    % 根据体系大小选择计算模式
    if N > 500
        % 大型体系: 使用并行计算加速
        fprintf('  体系较大(%d原子)，启用并行计算...\n', N);
        parfor k = 1:N_neighbor
            H_slice = zeros(N, N);  % 当前近邻的哈密顿量切片
            dx0 = disp_vec(k, 1);   % 当前近邻的x位移
            dy0 = disp_vec(k, 2);   % 当前近邻的y位移
            
            % 遍历所有原子对
            for i = 1:N
                for j = 1:N
                    % 计算原子间相对位置
                    Rx = x(i) - x(j) + dx0;
                    Ry = y(i) - y(j) + dy0;
                    Rz = z(i) - z(j);
                    R_sq = Rx^2 + Ry^2;  % 距离平方
                    R_ij = sqrt(R_sq + Rz^2);
                                        
                    % 跃迁条件：距离在合理范围内或不同层
                    if layer(i) == layer(j) && R_ij > 0.9*a_0 && R_ij < 1.1*a_0
                        H_slice(i, j) = H_slice(i, j) + t1;
                    elseif layer(i) ~= layer(j) && R_sq < 0.1
                        H_slice(i, j) = H_slice(i, j) + t2;
                    elseif i == j && k == 1 
                        if layer(i) == 1
                            H_slice(i, j) = H_slice(i, j) + E_elec;
                        elseif layer(i) == 2
                            H_slice(i, j) = H_slice(i, j) - E_elec;
                        end
                    end
                end
            end
            
            % 存储当前近邻的哈密顿量切片
            Hamiltonian(:, :, k) = H_slice;
        end
    else
        % 小型体系: 串行计算
        fprintf('  体系较小(%d原子)，使用串行计算...\n', N);
        for k = 1:N_neighbor
            dx0 = disp_vec(k, 1);  % 当前近邻的x位移
            dy0 = disp_vec(k, 2);  % 当前近邻的y位移
            
            % 遍历所有原子对
            for i = 1:N
                for j = 1:N
                    % 计算原子间相对位置
                    Rx = x(i) - x(j) + dx0;
                    Ry = y(i) - y(j) + dy0;
                    Rz = z(i) - z(j);
                    R_sq = Rx^2 + Ry^2;  % 距离平方
                    R_ij = sqrt(R_sq + Rz^2);
                    
                    % 跃迁条件：距离在合理范围内或不同层
                    if layer(i) == layer(j) && R_ij > 0.9*a_0 && R_ij < 1.1*a_0
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + t1;
                    elseif layer(i) ~= layer(j) && R_sq < 0.1
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + t2;
                    elseif i == j && k == 1 
                        if layer(i) == 1
                            Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + E_elec;
                        elseif layer(i) == 2
                            Hamiltonian(i, j, k) = Hamiltonian(i, j, k) - E_elec;
                        end
                    end
                end
            end
        end
    end
end

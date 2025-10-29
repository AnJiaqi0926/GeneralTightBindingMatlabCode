function [Hamiltonian,neighbor_cell, disp_vec] = BuildHamiltonianGraphene_SlaterKoster(a1, a2, atom_position, N, neighbor_cell, params)
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
    d_cut = (3 * a_lattice); % 截断长度=3倍石墨烯晶格常数
    d_cut_min = (0.01 * a_lattice); % 截断长度=3倍石墨烯晶格常数
    d_0 = 3.35;
    r0 = 0.184*a_lattice;
    

    % ================== S-K 参数 (GBN) ==================
    if isfield(params, 'Vpp_pi0')
        Vpp_pi0 = params.Vpp_pi0;
    else
        Vpp_pi0 = -2.7;
    end

    if isfield(params, 'Vpp_sigma0')
        Vpp_sigma0 = params.Vpp_sigma0;
    else
        Vpp_sigma0 = 0.48;
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
                    if R_sq < d_cut && (i ~= j || k ~= 1)
                        Vpp_pi = Vpp_pi0 * exp(-(R_ij - a_0)/r0);
                        Vpp_sigma = Vpp_sigma0 * exp(-(R_ij - d_0)/r0);
                        
                        % 计算Slater-Koster跃迁积分
                        t_val = Vpp_pi * (1 - (abs(Rz)/R_ij)^2) + Vpp_sigma * ((abs(Rz)/R_ij)^2);
                        H_slice(i, j) = H_slice(i, j) + t_val;
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
                    if R_sq < d_cut && R_sq > d_cut_min && (i ~= j || k ~= 1)
                        % 区分层内和层间跃迁
                        Vpp_pi = Vpp_pi0 * exp(-(R_ij - a_0)/r0);
                        Vpp_sigma = Vpp_sigma0 * exp(-(R_ij - d_0)/r0);
                        
                        % 计算Slater-Koster跃迁积分
                        t_val = Vpp_pi * (1 - (abs(Rz)/R_ij)^2) + Vpp_sigma * ((abs(Rz)/R_ij)^2);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + t_val;
                    end
                end
            end
        end
    end
end

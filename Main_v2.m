% ========================================================================
% 二维材料能带结构计算程序 (石墨烯示例)
% 功能：
%   1. 构建任意尺寸的石墨烯超胞
%   2. 计算哈密顿量矩阵
%   3. 支持两种能带计算模式：
%       - 全布里渊区模式
%       - 沿高对称线路径模式 (Γ-M-K-Γ)
%   4. 可视化能带结构和原子结构
% 
% 作者：An Jiaqi & DeepSeek
% 日期：2025年7月04日
% ========================================================================
% 接下来要做的：1，尝试读取LAMMPS输出的结构文件，用于计算BGBN等体系，2. 正确求解Displacement，修改SK GBN 哈密顿函数
% 3. 添加DOS calculat
% ========================================================================

clear all  % 清除所有变量
tic        % 开始计时

%% ========================== 初始设置(计算参数) ==========================
% ------------------------- 通用参数 -------------------------
nsum = 500;  % k点数量 (全布里渊区模式为nsum×nsum个点，高对称线模式为nsum个点)
nsum_chern = 50;  % 计算陈数的k点数量

% 态密度计算参数
delta_v = 1*1e-4;
omega = linspace(-3, 3, 1000);  % 能量范围 ω
d = 2;  % 假设是2D布里渊区
m = 500; % 扩展后的矩阵大小

% ----------------------- 计算模式设定 -----------------------
% 计算模式选择 ('full' 或 'high_sym')
calc_mode = 'high_sym';  % 'full' - 全布里渊区能带, 'high_sym' - 沿高对称线能带

calc_chern = false; % 是否计算陈数
calc_ribbon = false; % 是否计算条带结构
calc_DOS = true; % 是否计算态密度

calc_Structure = 'GBN_Moire'; % 计算Bulk结构 Single_Graphene Bilayer_Graphene GBN_Commensurate_Test GBN_Moire
calc_Structure_ribbon = 'zigzag_SingleGraphene';  % 计算ribbon结构 Armchair zigzag

calc_Hamiltonian = 'Slater_Koster_GBN'; % 计算Hamiltonian的方式 Slater_Koster_Graphene SimplestNearest_Graphene Test Slater_Koster_for_GBN

% ------------------------ 原胞构建参数 ----------------------
% 这些参数在石墨烯结构中未使用，保留作为通用接口

Nx = 5;  % 沿a1方向的超胞重复次数
Ny = 5;  % 沿a2方向的超胞重复次数
Ny_ribbon = 200; % 沿a2方向的ribbon超胞重复次数

plot_structure = true;  % 是否绘制原子结构图

% 定义9个最近邻原胞 (包括原胞自身)
neighbor_cell = [0, 0;   % 自身
                     1, 0;   % 正a1方向
                     0, 1;   % 正a2方向
                    -1, 1;   % 负a1正a2方向
                    -1, 0;   % 负a1方向
                     0,-1;   % 负a2方向
                     1,-1;   % 正a1负a2方向
                    -1,-1;   % 负a1负a2方向
                     1, 1];  % 正a1正a2方向

% 定义3个最近邻原胞 (包括原胞自身)
neighbor_cell_ribbon = [0, 0;   % 自身
                     1, 0;   % 正周期方向
                     -1, 0];  % 负周期方向

% ---------------------- Hamiltonian参数 ---------------------

%% ================== 结构构建函数映射 ==================
% 构建结构模式选择 ('Single_Graphene' 或 'GBN')
% 构建Hamiltonian模式选择 ('SimplestNearest', 'Slater_Koster', 'Slater_Koster for GBN')
structureFuncs = struct(...
    'Single_Graphene', @Structure_SingleGrapheneSupercell, ...
    'Bilayer_Graphene', @Structure_BilayerGrapheneSupercell, ...
    'GBN_Commensurate_Test', @Structure_GBNCommensuratecell, ...
    'GBN_Moire', @Structure_GBNMoirecell, ...
    'Armchair', @Structure_ArmchairGrapheneSupercell, ...
    'zigzag_BilayerGraphene',@Structure_zigzagGrapheneSupercell, ...
    'zigzag_SingleGraphene', @Structure_zigzagSingleGraphene);

structureFuncs2 = struct(...
    'ReadXYZ', @Structure_ReadXYZ, ...
    'BilayerG_on_BN_Build', @Structure_BGBNMoirecell);

hamiltonianFuncs = struct(...
    'Slater_Koster_Graphene', @BuildHamiltonianGraphene, ...
    'SimplestNearest_Graphene', @BuildHamiltonianGraphene_SimplestNearest, ...
    'GrapheneQAHE', @BuildHamiltonianGraphene_RashbaSOC_exchangefield_electricfield, ...
    'Slater_Koster_GBN', @BuildHamiltonianGBN);

%% ================== 构建几何结构 ==================
% 构建超胞
% 输出:
%   a1, a2 - 超胞基矢 (Å)
%   atom_position - 原子位置矩阵 [N×10]
%     9列分别对应: x, y, z, layer, sublattice, atom_type, spin, dx, dy, dz

% 检查是否支持该结构
if isfield(structureFuncs, calc_Structure)
    [a1, a2, atom_position] = structureFuncs.(calc_Structure)(Nx, Ny, plot_structure);
elseif isfield(structureFuncs2, calc_Structure)
    [a1, a2, atom_position] = structureFuncs2.(calc_Structure)(Nx, Ny, plot_structure);
else
    error('未知计算结构: %s', calc_Structure);
end

N = size(atom_position, 1);  % 体系中的原子总数
fprintf('系统总原子数: %d\n', N);

if calc_ribbon
    if isfield(structureFuncs, calc_Structure_ribbon)
        [a1_ribbon, a2_ribbon, atom_position_ribbon] = structureFuncs.(calc_Structure_ribbon)(1, Ny_ribbon, plot_structure);
    else
        error('未知计算结构: %s', calc_Structure_ribbon);
    end
    
    N_ribbon = size(atom_position_ribbon, 1);  % 体系中的原子总数
    fprintf('系统总原子数: %d\n', N_ribbon);
end

%% ================== 动量空间设置与k点生成 ==================
% 生成k点路径
% 输出:
%   b1, b2 - 倒格矢
%   kx_v, ky_v - k点坐标 (笛卡尔坐标)
%   kpath - 高对称线路径长度坐标 (仅在高对称线模式下有效)
%   high_sym_points - 高对称点分数坐标
[b1, b2, kx_v, ky_v, kpath, high_sym_points] = ChooseKpoint(a1, a2, nsum, calc_mode);

% 确定总k点数
nsum1 = size(kx_v, 2);  % 实际计算的k点数量
fprintf('计算模式: %s, k点数: %d\n', calc_mode, nsum1);

% 初始化能带存储矩阵
Energy_v = zeros(N, nsum1);  % 存储每个k点的能量本征值

%% ================== Ribbon的动量空间设置与k点生成 ==================
% 生成k点路径
% 输出:
%   b1, b2 - 倒格矢
%   kx_v, ky_v - k点坐标 (笛卡尔坐标)
%   kpath - 高对称线路径长度坐标 (仅在高对称线模式下有效)
%   high_sym_points - 高对称点分数坐标
if calc_ribbon
    [b1_ribbon, b2_ribbon, kx_v_ribbon] = ChooseKpoint_Ribbon(a1_ribbon, a2_ribbon, nsum);
    kpath_ribbon = kx_v_ribbon;
    
    % 确定总k点数
    nsum1_ribbon = size(kx_v_ribbon, 2);  % 实际计算的k点数量
    fprintf('计算模式: %s, k点数: %d\n', calc_mode, nsum1_ribbon);
    
    % 初始化能带存储矩阵
    Energy_v_ribbon = zeros(N_ribbon, nsum1_ribbon);  % 存储每个k点的能量本征值
end

%% ================== 构建哈密顿量矩阵 ==================
% 构建紧束缚模型哈密顿量
% 输出:
%   Hamiltonian - [N×N×N_neighbor] 哈密顿量张量
%   neighbor_cell - 近邻原胞列表
%   disp_vec - 位移矢量
if isfield(hamiltonianFuncs, calc_Hamiltonian)
    [Hamiltonian, neighbor_cell, disp_vec] = hamiltonianFuncs.(calc_Hamiltonian)(a1, a2, atom_position, N, neighbor_cell);
else
    error('不支持的哈密顿量结构: %s', calc_Hamiltonian);
end

N_neighbor = size(neighbor_cell, 1);  % 近邻原胞数量
fprintf('哈密顿量构建完成, 近邻原胞数: %d\n', N_neighbor);

if calc_ribbon
    if isfield(structureFuncs, calc_Structure_ribbon)
        [Hamiltonian_ribbon, neighbor_cell_ribbon, disp_vec_ribbon] = hamiltonianFuncs.(calc_Hamiltonian)(a1_ribbon, a2_ribbon, atom_position_ribbon, N_ribbon, neighbor_cell_ribbon);
    else
        error('不支持的哈密顿量结构: %s', calc_Hamiltonian);
    end
    
    N_neighbor_ribbon = size(neighbor_cell_ribbon, 1);  % 近邻原胞数量
    fprintf('哈密顿量构建完成, 近邻原胞数: %d\n', N_neighbor_ribbon);
end

%% ================== 能带结构计算 ==================
fprintf('开始能带计算...\n');
for i = 1:nsum1
    kx = kx_v(i);  % 当前k点的x分量
    ky = ky_v(i);  % 当前k点的y分量
    
    % 初始化k点相关哈密顿量
    HamiltonianK = zeros(N, N);
    
    % 对每个近邻原胞求和
    for k = 1:N_neighbor
        % 计算相位因子: e^{-i k·R}
        phase = exp(-1i*(kx*disp_vec(k,1) + ky*disp_vec(k,2)));
        
        % 添加近邻贡献
        HamiltonianK = HamiltonianK + Hamiltonian(:,:,k) * phase;
    end
    
    % 确保哈密顿量厄米性
    HamiltonianK = (HamiltonianK + HamiltonianK')/2;
    
    % 对角化哈密顿量，得到能量本征值
    Energy = eig(HamiltonianK);
    
    % 存储能量本征值
    Energy_v(:, i) = Energy;
    
    % 显示进度
    if mod(i, 100) == 0
        fprintf('  已完成 %d/%d k点\n', i, nsum1);
    end
end
fprintf('能带计算完成!\n');

%% ==================Ribbon 能带结构计算 ==================
if calc_ribbon
    fprintf('开始Ribbon能带计算...\n');
    for i = 1:nsum1_ribbon
        kx = kx_v_ribbon(i);  % 当前k点的x分量
        
        % 初始化k点相关哈密顿量
        HamiltonianK_ribbon = zeros(N_ribbon, N_ribbon);
        
        % 对每个近邻原胞求和
        for k = 1:N_neighbor_ribbon
            % 计算相位因子: e^{-i k·R}
            phase = exp(-1i*(kx*disp_vec_ribbon(k,1) ));
            
            % 添加近邻贡献
            HamiltonianK_ribbon = HamiltonianK_ribbon + Hamiltonian_ribbon(:,:,k) * phase;
        end
        
        % 确保哈密顿量厄米性
        HamiltonianK_ribbon = (HamiltonianK_ribbon + HamiltonianK_ribbon')/2;
        
        % 对角化哈密顿量，得到能量本征值
        Energy_ribbon = eig(HamiltonianK_ribbon);
        
        % 存储能量本征值
        Energy_v_ribbon(:, i) = Energy_ribbon;
        
        % 显示进度
        if mod(i, 100) == 0
            fprintf('  已完成 %d/%d k点\n', i, nsum1_ribbon);
        end
    end
    fprintf('能带计算完成!\n');
end

%% ================== 陈数计算 态密度计算  ==================
if calc_chern || calc_DOS
    fprintf('开始陈数与态密度计算...\n');
    
    % 设置陈数计算参数
    Nk = nsum_chern; % k点网格大小
    target_band = 1:round(N/2); % 计算费米面以下的能带
    band = target_band; % 目标能带索引
    
    % 准备k点网格
    kx_grid = linspace(-b1(1), b1(1), Nk);
    ky_grid = linspace(-b2(2), b2(2), Nk);
    dkx = kx_grid(2) - kx_grid(1); % kx微小偏移量
    dky = ky_grid(2) - ky_grid(1); % ky微小偏移量
    [KX, KY] = meshgrid(kx_grid, ky_grid);
    kx_v = KX(:)';
    ky_v = KY(:)';
    
    % 调用陈数计算函数

    [Omega, Chern_number, Chern_v1, Chern_v2, Ene_v] = ...
        Chern_Calculation(band, Hamiltonian, disp_vec, Nk, N, kx_grid, ky_grid, dkx, dky, b1, b2);

    if calc_chern
        % 输出结果
        fprintf('==============================================\n');
        fprintf('陈数计算结果: %.6f\n', Chern_number);
        fprintf('谷1陈数: %.6f, 谷2陈数: %.6f\n', Chern_v1, Chern_v2);
        fprintf('==============================================\n');
        
        % 绘制Berry曲率分布
        figure('Name', 'Berry曲率分布', 'Position', [100, 100, 800, 600]);
        imagesc(reshape(Omega, Nk, Nk));
        colorbar;
        if size(target_band,2) < 1
            title(sprintf('Berry曲率分布 (第%d能带)', target_band));
        elseif size(target_band,2) == round(N/2)
            title(sprintf('Berry曲率分布 (导带)'));
        end
        xlabel('k_x');
        ylabel('k_y');
        axis equal tight;
    end

    if calc_DOS
        Hsize = N;
        nsum1 = Nk;
                
        rho = zeros(size(omega));  % 初始化态密度
        for i = 1:Hsize
            % 创建原始网格
            x_old = linspace(1, nsum1, nsum1);  % 从1到n的n个点
            y_old = linspace(1, nsum1, nsum1);  % 从1到n的n个点
            
            % 创建目标网格
            x_new = linspace(1, nsum1, m);  % 从1到n的m个点
            y_new = linspace(1, nsum1, m);  % 从1到n的m个点

            Ene_v_band = Ene_v(:,:,i);
            
            % 使用二维插值
            [Ene_new] = interp2(x_old, y_old, Ene_v_band, x_new', y_new, 'cubic');

            for i_rho = 1:length(omega)
                % 计算 δ(ω - ϵμ(k))，这里用 Dirac 近似 (即通过非常小的值来模拟 δ 函数)
                delta = @(omega_k) (abs(omega_k - omega(i_rho)) < delta_v);  % 假设 δ 函数近似为一个小的范围
                
                % 对布里渊区进行积分（数值积分）
                kx = linspace(-pi, pi, m);  % x方向的k点
                ky = linspace(-pi, pi, m);  % y方向的k点
                
                [KX, KY] = meshgrid(kx, ky);  % 创建k空间网格
                % 从能带数据中提取对应的能量
                epsilon_k = Ene_new;  % 取出第mu条能带的能量数据
        
                % 计算 δ(ω - ϵμ(k))，并对k空间进行积分
                integrand = delta(epsilon_k);  % δ函数的值
                
                % 计算积分，这里使用了矩形积分近似
                rho(i_rho) = rho(i_rho) + sum(integrand(:)) * (2 * pi)^(-d) / m^2;  % 累加积分结果
            end
            if mod(i, 10) == 0
                fprintf('  已完成 %d/%d k点\n', i, Hsize);
            end
        end

        % 绘制结果
        figure;
        plot(omega,rho);
        xlim([-3,3])
        xlabel('ω');
        ylabel('ρ(ω)');
        title('态密度 ρ(ω) 的计算');
        grid on;
        
        omega = omega';
        rho = rho';
        
        save('Ene_0.1GPa.dat', 'omega', '-ascii');
        save('DOS_0.1GPa.dat', 'rho', '-ascii');
    end
        
end

%% ================== 结果可视化 ==================
if strcmp(calc_mode, 'high_sym')
    % 保存能带数据
    plot_dat = [kpath(:), Energy_v'];
    save plot_band_high_sym.dat plot_dat -ascii;
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 1;  % 能带数量
    
    kpath(:) = plot_dat(:,1);  % kpath坐标

    % ========== 高对称线模式绘图 ========== 
    fprintf('绘制高对称线能带结构...\n');
    figure('Name', '高对称线能带结构', 'Position', [100, 100, 800, 600]);

    if Hsize > 40
        % 绘制费米面附近的40条能带
        start_band = max(1, round(Hsize/2) - 20);
        end_band = min(Hsize, round(Hsize/2) + 19);
        num_bands = 40;
        fprintf('绘制费米面附近的%d条能带 (从第%d条到第%d条)\n', num_bands, start_band, end_band);
    else
        % 能带总数小于等于40，绘制全部能带
        start_band = 1;
        end_band = Hsize;
        num_bands = Hsize;
        fprintf('能带总数较少(%d条)，绘制全部能带\n', Hsize);
    end
    
    % 绘制所有能带
    for i = start_band:end_band
        plot(kpath, Energy_v(i, :), 'b-', 'LineWidth', 0.5);
        hold on;
    end
    
    % 标记高对称点位置
    sym_points = unique(kpath([1, round(nsum1/3)+1, round(2*nsum1/3)+1, end]));
    for i = 1:length(sym_points)
        xline(sym_points(i), '--k', 'LineWidth', 1.5);
    end
    
    % 添加高对称点标签
    labels = {'\Gamma', 'M', 'K', '\Gamma'};  % 高对称点标签
    y_min = min(Energy_v(:));  % 能量最小值
    y_max = max(Energy_v(:));  % 能量最大值
    y_range = y_max - y_min;   % 能量范围
    
    % 在图中标注高对称点
    text(sym_points(1), y_min - 0.05*y_range, labels{1}, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    text(sym_points(2), y_min - 0.05*y_range, labels{2}, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    text(sym_points(3), y_min - 0.05*y_range, labels{3}, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    text(sym_points(4), y_min - 0.05*y_range, labels{4}, ...
        'FontSize', 14, 'HorizontalAlignment', 'center');
    
    % 设置图形属性
    title(sprintf('%dx%d 能带结构 (沿 \\Gamma-M-K-\\Gamma 路径)', Nx, Ny));
    xlabel('K path');
    ylabel('Energy (eV)');
    grid on;
    set(gca, 'FontSize', 12);
    set(gca, 'Box', 'on');
    
elseif strcmp(calc_mode, 'full')
    % ========== 全布里渊区模式绘图 ========== 
    fprintf('绘制全布里渊区能带结构...\n');
    
    % 保存能带数据
    plot_dat = [kx_v(:), ky_v(:), Energy_v'];
    save plot_band_full.dat plot_dat -ascii;
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 2;  % 能带数量
    
    kx = plot_dat(:,1);  % kx坐标
    ky = plot_dat(:,2);  % ky坐标
    
    % 创建新图形
    figure('Name', 'Band Structure', 'Position', [100, 100, 1000, 800]);

    if Hsize > 40
        % 绘制费米面附近的40条能带
        start_band = max(1, round(Hsize/2) - 20);
        end_band = min(Hsize, round(Hsize/2) + 19);
        num_bands = 40;
        fprintf('绘制费米面附近的%d条能带 (从第%d条到第%d条)\n', num_bands, start_band, end_band);
    else
        % 能带总数小于等于40，绘制全部能带
        start_band = 1;
        end_band = Hsize;
        num_bands = Hsize;
        fprintf('能带总数较少(%d条)，绘制全部能带\n', Hsize);
    end
    
    % 绘制每个能带的3D曲面
    for band_idx = start_band:end_band  % 最多绘制费米面附近40能带
        % 创建网格数据
        kx_v_grid = linspace(min(kx), max(kx), nsum);
        ky_v_grid = linspace(min(ky), max(ky), nsum);
        Ene_v = zeros(nsum, nsum);  % 能量矩阵
        
        % 填充能量矩阵
        index = 0;
        for ii = 1:nsum
            for jj = 1:nsum
                index = index + 1;
                Ene_v(ii, jj) = plot_dat(index, band_idx + 2);
            end
        end
        
        % 绘制3D曲面
        surf(kx_v_grid, ky_v_grid, Ene_v);
        hold on;
    end
    
    % 设置图形属性
    title(sprintf('%dx%d 石墨烯超胞全布里渊区能带结构', Nx, Ny));
    xlabel('k_x');
    ylabel('k_y');
    zlabel('Energy (eV)');
    grid on;
    colormap(jet);  % 使用彩色映射
    colorbar;       % 添加颜色条
    view(45, 30);   % 设置视角
end

%% ================== Ribbon 结果可视化 ==================
if calc_ribbon
    % 保存能带数据
    plot_dat = [kpath_ribbon(:), Energy_v_ribbon'];
    save plot_band_ribbon.dat plot_dat -ascii;
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 1;  % 能带数量
    
    kpath(:) = plot_dat(:,1);  % kpath坐标

%% ========== Ribbon 能带绘图 ========== 

    fprintf('绘制高对称线能带结构...\n');
    figure('Name', 'Ribbon 能带结构', 'Position', [100, 100, 800, 600]);
    
    % 能带总数小于等于40，绘制全部能带
    start_band = 1;
    end_band = Hsize;
    num_bands = Hsize;
    fprintf('能带总数较少(%d条)，绘制全部能带\n', Hsize);
    
    % 绘制所有能带
    for i = start_band:end_band
        plot(kpath, Energy_v_ribbon(i, :), 'b-', 'LineWidth', 0.5);
        hold on;
    end
    
    % 设置图形属性
    title(sprintf('%dx%d 能带结构 (沿 \\Gamma-M-K-\\Gamma 路径)', Nx, Ny));
    xlabel('K path');
    ylabel('Energy (eV)');
    grid on;
    set(gca, 'FontSize', 12);
    set(gca, 'Box', 'on');
end


toc  % 显示总计算时间
fprintf('程序执行完成!\n');
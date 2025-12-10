 % ========================================================================
% 二维材料能带结构计算程序 (石墨烯示例)
% 功能：
%   1. 构建任意尺寸的二维材料超胞
%   2. 计算紧束缚模型哈密顿量矩阵
%   3. 支持两种能带计算模式：
%       - 全布里渊区模式
%       - 沿高对称线路径模式 (Γ-M-K-Γ)
%   4. 可视化能带结构和原子结构
% 
% 作者：An Jiaqi & DeepSeek
% 日期：2025年7月04日
% ========================================================================
% 正在优化的事： 1. 
%               2. 
% ========================================================================
% 待要优化的事：1. 
% ========================================================================

clear all  % 清除所有变量
tic        % 开始计时

%% ================== 结构构建函数映射 ==================
structureFuncs = struct(...
    'GBN', @Structure_GBN, ...
    'GBN_ReadXYZ', @Structure_GBN_ReadXYZ, ...
    'Multilayer_Graphene', @Structure_MultilayerGraphene, ...
    'zigzag_Multilayer_Graphene', @Structure_zigzagMultilayerGraphene);

hamiltonianFuncs = struct(...
    'Slater_Koster_GBN', @BuildHamiltonianGBNSlaterKoster, ...
    'GBN_FullTB', @BuildHamiltonianGBN, ...
    'SimpleNearest', @BuildHamiltonianSimpleNearest);

%% ========================== 初始设置(计算参数) ==========================
% ----------------------- 计算模式设定 -----------------------
% 计算模式选择 ('full' 或 'high_sym')
calc_Structure = 'GBN_ReadXYZ'; % 计算Bulk结构
calc_Structure_ribbon = 'zigzag_Multilayer_Graphene';  % 计算ribbon结构

calc_Hamiltonian = 'Slater_Koster_GBN'; % 计算Hamiltonian的方式

calc_mode = 'high_sym';  % 'full' - 全布里渊区能带, 'high_sym' - 沿高对称线能带
calc_ribbon = false; % 是否计算条带结构
calc_chern = false; % 是否计算陈数
calc_DOS = false; % 是否计算态密度

% ------------------------ 原胞构建参数 ----------------------
Nx = 1;  % 沿a1方向的超胞重复次数
Ny = 1;  % 沿a2方向的超胞重复次数
Nx_ribbon = 1; %沿a1方向的ribbon超胞重复次数
Ny_ribbon = 100; % 沿a2方向的ribbon超胞重复次数

plot_structure = false;  % 是否绘制原子结构图

% 在Main_TB.m中根据体系类型调整neighbor_cell
% 原有的9近邻定义
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

% 定义ribbon的3个最近邻原胞 (包括原胞自身)
neighbor_cell_ribbon = [0, 0;   % 自身
                        1, 0;   % 正周期方向
                       -1, 0];  % 负周期方向

structure_params = GetStructureParams(calc_Structure);
structure_params_ribbon = GetStructureParams(calc_Structure_ribbon);

% ---------------------- Hamiltonian参数 ---------------------
hamiltonian_params = GetHamiltonianParams(calc_Hamiltonian);

% ------------------------- 态密度计算参数 -------------------------
delta_v = 1*1e-2; % 态密度使用的delta函数展宽
d = 1;  % 布里渊区面积
numInterpPoints = 2000; % 插值法绘制3D能带

% ------------------------- 其他通用参数 -------------------------
nsum = 1001;  % 绘制能带的k点数量 (全布里渊区模式为nsum×nsum个点，高对称线模式为nsum个点)

% custom_path = [0, 0;       % Γ点
%                0, 1/2;     % M点  
%                1/3, 2/3;   % K点
%                0, 0];      % Γ点 (闭合路径)

custom_path = [0, 0;       % Γ点
               1/3, 2/3;     % M点  
               2/3, 4/3;   % K点
               1, 2];      % Γ点 (闭合路径)

% ---------------------- AAH模型特定参数---------------
calc_time_evolution = false;  % 是否计算含时演化
time_range = [0, 10];        % 时间范围
dt = 0.01;                   % 时间步长

% ---------------------- GBN模型特定参数---------------
if strcmp(calc_Structure, 'GBN_ReadXYZ')
    [value,Name] = textread('Input.in','%f%s');
    
    layer2 = value(2);
    layer1 = value(3);
    
    aG = value(4);
    aBN = value(5);
    
    theta = value(6);
    phi = value(7);

    structure_params.layer2 = value(2);
    structure_params.layer1 = value(3);
end

%% ================== 构建几何结构 ==================
% 构建超胞
% 输出:
%   a1, a2 - 超胞基矢 (Å)n 
%   atom_position - 原子位置矩阵 [N×10]
%     9列分别对应: x, y, z, layer, sublattice, atom_type, spin, dx, dy, dz

% 检查是否支持该结构
if isfield(structureFuncs, calc_Structure)
    if strcmp(calc_Structure, 'GBN_ReadXYZ')
        [a1, a2, atom_position,atom_position_Init] = structureFuncs.(calc_Structure)(Nx, Ny, plot_structure,structure_params);
    else
        [a1, a2, atom_position] = structureFuncs.(calc_Structure)(Nx, Ny, plot_structure,structure_params);
    end
else
    error('未知计算结构: %s', calc_Structure);
end

N = size(atom_position, 1);  % 体系中的原子总数
fprintf('系统总原子数: %d\n', N);

if calc_ribbon
    if isfield(structureFuncs, calc_Structure_ribbon)
        [a1_ribbon, a2_ribbon, atom_position_ribbon] = structureFuncs.(calc_Structure_ribbon)(Nx_ribbon, Ny_ribbon, plot_structure,structure_params_ribbon);
    else
        error('未知Ribbon计算结构: %s', calc_Structure_ribbon);
    end
    
    N_ribbon = size(atom_position_ribbon, 1);  % 体系中的原子总数
    fprintf('Ribbon结构系统总原子数: %d\n', N_ribbon); 
end

%% ==================Bulk 动量空间设置与k点生成 ==================
% 生成k点路径
% 输出:
%   b1, b2 - 倒格矢
%   kx_v, ky_v - k点坐标 (笛卡尔坐标)
%   kpath - 高对称线路径长度坐标 (仅在高对称线模式下有效)
%   high_sym_points - 高对称点分数坐标
[b1, b2, kx_v, ky_v, kpath, high_sym_points] = ChooseKpoint(a1, a2, nsum, calc_mode,custom_path);

% 确定总k点数
nsum1 = size(kx_v, 2);  % 实际计算的k点数量
fprintf('计算模式: %s, k点数: %d\n', calc_mode, nsum1);

% 初始化能带存储矩阵
Energy_v = zeros(N, nsum1);  % 存储每个k点的能量本征值

% ================== Ribbon的动量空间设置与k点生成 ==================
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
    fprintf('Ribbon 能带计算的k点数: %d\n', nsum1_ribbon);
    
    % 初始化能带存储矩阵
    Energy_v_ribbon = zeros(N_ribbon, nsum1_ribbon);  % 存储每个k点的能量本征值
end

%% ================== 构建哈密顿量矩阵 ==================
% 构建紧束缚模型哈密顿量
% 输出:
%   Hamiltonian
%   neighbor_cell - 近邻原胞列表
%   disp_vec - 位移矢量
% 获取参数

if isfield(hamiltonianFuncs, calc_Hamiltonian)
    if strcmp(calc_Hamiltonian, 'GBN_FullTB')
        [Hamiltonian, neighbor_cell, disp_vec] = hamiltonianFuncs.(calc_Hamiltonian)(a1, a2, atom_position,atom_position, N, neighbor_cell, hamiltonian_params);
    else
        [Hamiltonian, neighbor_cell, disp_vec] = hamiltonianFuncs.(calc_Hamiltonian)(a1, a2, atom_position, N, neighbor_cell, hamiltonian_params);
    end
else
    error('不支持的哈密顿量结构: %s', calc_Hamiltonian);
end

N_neighbor = size(neighbor_cell, 1);  % 近邻原胞数量
%fprintf('Bulk哈密顿量构建完成, 近邻原胞数: %d\n', N_neighbor);

if calc_ribbon
    if isfield(structureFuncs, calc_Structure_ribbon)
        [Hamiltonian_ribbon, neighbor_cell_ribbon, disp_vec_ribbon] = hamiltonianFuncs.(calc_Hamiltonian)(a1_ribbon, a2_ribbon, atom_position_ribbon, N_ribbon, neighbor_cell_ribbon, hamiltonian_params);
    else
        error('不支持的哈密顿量结构: %s', calc_Hamiltonian);
    end
    
    N_neighbor_ribbon = size(neighbor_cell_ribbon, 1);  % 近邻原胞数量
    fprintf('Ribbon哈密顿量构建完成, 近邻原胞数: %d\n', N_neighbor_ribbon);
end

%% ================== Bulk 能带结构计算 ==================
fprintf('开始Bulk能带计算...\n');
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
    % if mod(i, 100) == 0
    %    fprintf('  已完成 %d/%d k点\n', i, nsum1);
    % end
    fprintf('  已完成 %d/%d k点\n', i, nsum1);
end
fprintf('Bulk能带计算完成!\n');

plot_dat = [kpath(:), Energy_v'];
filename = "Band_" + calc_Structure + "_" + calc_Hamiltonian + ".dat";
save(filename, "plot_dat", "-ascii");

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
    fprintf('Ribbon能带计算完成!\n');
end

%% ================== 结果可视化 ==================
if strcmp(calc_mode, 'high_sym')
    % 保存能带数据
    plot_dat = [kpath(:), Energy_v'];
    save("Band_"+calc_Structure+"_high_sym.dat","plot_dat","-ascii");
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 1;  % 能带数量
    
    kpath(:) = plot_dat(:,1);  % kpath坐标

    % ========== 高对称线模式绘图 ========== 
    fprintf('绘制高对称线能带结构...\n');
    fig_high_sym = figure('Name', '高对称线能带结构', 'Position', [100, 100, 800, 600]);

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
    % sym_points = unique(kpath([1, round(nsum1/3)+1, round(2*nsum1/3)+1, end]));
    % 准确计算高对称点在kpath中的位置
    sym_points = zeros(1, size(high_sym_points, 1));
    cumulative_length = 0;
    n_segments = size(high_sym_points, 1) - 1;
    for seg_idx = 1:n_segments
        start_point = high_sym_points(seg_idx, :);
        end_point = high_sym_points(seg_idx + 1, :);
        
        % 转换为笛卡尔坐标计算距离
        start_cart = start_point(1) * b1 + start_point(2) * b2;
        end_cart = end_point(1) * b1 + end_point(2) * b2;
        
        segment_length = norm(end_cart - start_cart);
        sym_points(seg_idx + 1) = cumulative_length + segment_length;
        cumulative_length = cumulative_length + segment_length;
    end

    for i = 1:length(sym_points)
        xline(sym_points(i), '--k', 'LineWidth', 1.5);
    end
    
    % 添加高对称点标签
    labels = {'\Gamma', 'K', 'K\prime', '\Gamma'};  % 高对称点标签
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

    savefig(fig_high_sym, "band_"+calc_Structure+"_high_sym.fig"); % 添加保存指令
    fprintf('高对称线能带结构已保存为 band_high_sym.fig\n');
    
elseif strcmp(calc_mode, 'full')
    % ========== 全布里渊区模式绘图 ========== 
    fprintf('绘制全布里渊区能带结构...\n');
    
    % 保存能带数据
    plot_dat = [kx_v(:), ky_v(:), Energy_v'];
    save("Band_"+calc_Structure+"_full.dat","plot_dat","-ascii");
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 2;  % 能带数量
    
    kx = plot_dat(:,1);  % kx坐标
    ky = plot_dat(:,2);  % ky坐标
    
    % 创建新图形
    fig_full = figure('Name', 'Band Structure', 'Position', [100, 100, 800, 600]);

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

    % 将k点网格化
    kx_grid = reshape(kx, nsum, nsum);
    ky_grid = reshape(ky, nsum, nsum);
    
    % 绘制每个能带的3D曲面
    for band_idx = start_band:end_band
        % 提取当前能带的能量值并网格化
        Ene_band = plot_dat(:, band_idx + 2);
        Ene_grid = reshape(Ene_band, nsum, nsum);
        
        % 绘制3D曲面
        surf(kx_grid, ky_grid, Ene_grid);
        hold on;
    end
    
    % 设置图形属性
    title(sprintf('%dx%d 石墨烯超胞全布里渊区能带结构', Nx, Ny));
    xlabel('k_x');
    ylabel('k_y');
    zlabel('Energy (eV)');
    grid off;
    axis equal;
    colormap(jet);  % 使用彩色映射
    colorbar;       % 添加颜色条
    view(45, 30);   % 设置视角

    % 保存为FIG文件
    savefig(fig_full, "band_"+calc_Structure+"_full.fig"); % 添加保存指令
    fprintf('全布里渊区能带结构已保存为 band_full.fig\n');
end

%% ================== Ribbon 结果可视化 ==================
if calc_ribbon
    % 保存能带数据
    plot_dat = [kpath_ribbon(:), Energy_v_ribbon'];
    save("Band_"+calc_Structure+"_ribbon.dat","plot_dat","-ascii");
    
    % 准备绘图数据
    size1 = size(plot_dat);
    Hsize = size1(2) - 1;  % 能带数量
    
    kpath_ribbon(:) = plot_dat(:,1);  % kpath坐标

%% ========== Ribbon 能带绘图 ========== 

    fprintf('绘制高对称线能带结构...\n');
    fig_ribbon = figure('Name', 'Ribbon 能带结构', 'Position', [100, 100, 800, 600]);
    
    % 能带总数小于等于40，绘制全部能带
    start_band = 1;
    end_band = Hsize;
    num_bands = Hsize;
    fprintf('能带总数较少(%d条)，绘制全部能带\n', Hsize);
    
    % 绘制所有能带
    for i = start_band:end_band
        plot(kpath_ribbon, Energy_v_ribbon(i, :), 'b-', 'LineWidth', 0.5);
        hold on;
    end
    
    % 设置图形属性
    title(sprintf('%dx%d 能带结构 (沿 \\Gamma-M-K-\\Gamma 路径)', Nx, Ny));
    xlabel('K path');
    ylabel('Energy (eV)');
    grid on;
    set(gca, 'FontSize', 12);
    set(gca, 'Box', 'on');

    % 保存为FIG文件
    savefig(fig_ribbon, "band_"+calc_Structure+"_ribbon.fig"); % 添加保存指令
    fprintf('Ribbon能带结构已保存为 band_ribbon.fig\n');
end

%% ================== 态密度计算  ==================
if calc_DOS 
    fprintf('开始态密度计算...\n');
    bar = waitbar(0,'读取数据中...');
    if strcmp(calc_mode, 'full')
        Hsize = size(Energy_v);
        Hsize = Hsize(1);
        plot_dat = [kx_v(:), ky_v(:), Energy_v'];
    else
        fprintf('需要重新计算全布里渊区能带...\n');
        [b1, b2, kx_v, ky_v, kpath, high_sym_points] = ChooseKpoint(a1, a2, nsum, 'full',custom_path);
        % 确定总k点数
        nsum1 = size(kx_v, 2);  % 实际计算的k点数量
        fprintf('计算模式: %s, k点数: %d\n', calc_mode, nsum1);
        
        % 初始化能带存储矩阵
        Energy_v = zeros(N, nsum1);  % 存储每个k点的能量本征值

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
            str=['计算中...',num2str(round(100*i/nsum1)),'%'];
            waitbar(i/nsum1,bar,str) 
        end
        fprintf('态密度全布里渊区能带计算完成!\n');
        Hsize = size(Energy_v);
        Hsize = Hsize(1);
    end
    close(bar) 


    x_old = linspace(-0.5, 0.5, nsum);
    y_old = linspace(-0.5, 0.5, nsum);

    % 创建目标网格
    x_new = linspace(-0.5, 0.5, numInterpPoints);
    y_new = linspace(-0.5, 0.5, numInterpPoints);

    omega = linspace(min(Energy_v(:))-abs(0.1*min(Energy_v(:))), max(Energy_v(:))+abs(0.1*max(Energy_v(:))), numInterpPoints);  % 态密度计算能量范围 ω
    rho = zeros(1,numInterpPoints);  % 初始化态密度

    for i = 1:Hsize
        Ene_v_band = reshape(Energy_v(i,:),nsum,nsum);

        % 使用二维插值
        [Ene_new] = interp2(x_old, y_old, Ene_v_band, x_new', y_new, 'cubic');

        for i_rho = 1:numInterpPoints
            % 计算 δ(ω - ϵμ(k))，这里用 Dirac 近似 (即通过非常小的值来模拟 δ 函数)
            delta = @(omega_k) (abs(omega_k - omega(i_rho)) < delta_v);  % 假设 δ 函数近似为一个小的范围

            % 对布里渊区进行积分（数值积分）
            [KX, KY] = meshgrid(x_new,y_new);  % 创建k空间网格
            % 从能带数据中提取对应的能量
            epsilon_k = Ene_new;  % 取出第mu条能带的能量数据

            % 计算 δ(ω - ϵμ(k))，并对k空间进行积分
            integrand = delta(epsilon_k);  % δ函数的值

            % 计算积分，这里使用了矩形积分近似
            rho(i_rho) = rho(i_rho) + sum(integrand(:)) * (2 * pi)^(-d) / numInterpPoints^2;  % 累加积分结果
        end
        fprintf('  已完成 %d/%d 能量计算\n', i, Hsize);
    end

    % 绘制结果
    figure('Name', '态密度分布', 'Position', [100, 100, 800, 600]);
    plot(omega, rho, 'LineWidth', 2);
    xlim([min(omega), max(omega)]);
    xlabel('ω (eV)');
    ylabel('ρ(ω)');
    title('态密度 ρ(ω)');
    grid on;
    
    % 保存数据
    omega_col = omega(:);
    rho_col = rho(:);
    DOS = [omega_col, rho_col];
    save("DOS_"+calc_Structure+".dat", 'DOS', '-ascii');
    fprintf("态密度数据已保存到 DOS_"+calc_Structure+".dat\n");
    
    % 显示统计信息
    fprintf('态密度计算完成:\n');
    fprintf('  频率范围: %.3f 到 %.3f eV\n', min(omega), max(omega));
    fprintf('  最大态密度: %.6f\n', max(rho));
    fprintf('  使用的k点数: %d\n', numInterpPoints);
end

%% ================== 陈数计算  ==================
if calc_chern
    fprintf('开始陈数计算...\n');
    
    % 设置陈数计算参数
    Nk = nsum; % k点网格大小
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
        Chern_Calculation(band, Hamiltonian, disp_vec, Nk, N, kx_grid, ky_grid, dkx, dky, b1, b2,calc_Structure,calc_Hamiltonian);

    kx_grid_new = linspace(-b1(1), b1(1), numInterpPoints);
    ky_grid_new = linspace(-b2(2), b2(2), numInterpPoints);
    [Omega_new] = interp2(kx_grid, ky_grid, Omega, kx_grid_new', ky_grid_new, 'cubic');

    % 输出结果
    fprintf('==============================================\n');
    fprintf('陈数计算结果: %.6f\n', Chern_number);
    fprintf('谷1陈数: %.6f, 谷2陈数: %.6f\n', Chern_v1, Chern_v2);
    fprintf('==============================================\n');
    
    % 绘制Berry曲率分布
    figure('Name', 'Berry曲率分布', 'Position', [100, 100, 800, 600]);
    imagesc(reshape(Omega_new, numInterpPoints, numInterpPoints));
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

%% ================== 计算完成  ==================
toc  % 显示总计算时间
fprintf('程序执行完成!\n');

%% =================运算所需的子程序===================
function params = GetStructureParams(calc_Structure)
    % 根据哈密顿量类型和结构类型返回相应参数
    
    params = struct();
    
    switch calc_Structure
        case 'Multilayer_Graphene'
            params.ifspin = false;  % 是否考虑自旋
            params.bulkcellMode = "9cells"; % 近邻原胞模式
            params.a = 2.46;       % 石墨烯晶格常数 (Å)
            params.total_layer = 4;  % 结构的总层数
            params.StackedMode = 'ABC'; % 堆叠模式

        case 'zigzag_Multilayer_Graphene'
            params.ifspin = false;  % 是否考虑自旋
            params.bulkcellMode = "9cells"; % 近邻原胞模式
            params.a = 2.46;       % 石墨烯晶格常数 (Å)
            params.total_layer = 4;  % 结构的总层数
            params.StackedMode = 'ABC'; % 堆叠模式

        case 'GBN_ReadXYZ'
            params.Flag = 'ReadXYZ';
            params.layer1 = 1.0;
            params.layer2 = 1.0;
    end
end

function params = GetHamiltonianParams(calc_Hamiltonian)
    % 根据哈密顿量类型和结构类型返回相应参数
    
    params = struct();
    
    switch calc_Hamiltonian            
        case 'SimpleNearest'
            params.t1 = -2.7;
            params.t2 = 0.48;
            params.E_elec = 0.02;
            params.a = 2.46;
            % ... 其他参数

        case 'GBN_FullTB'
            [value,Name] = textread('Input.in','%f%s');

            params.layer2 = value(2);
            params.layer1 = value(3);
            
            params.aG = value(4);
            params.aBN = value(5);
            
            params.theta = value(6);
            params.phi = value(7);

            params.eps = params.aG/params.aBN-1;
            params.eps2 = params.aBN/params.aG-1;

            params.CAB = 1.987/1000;
            params.aG = 3.5*pi/180;
            params.CAB_BN = 4.418/1000;
            params.phiAB_BN = 26.10*pi/180;
            params.phase_mode = 'Nofull';

            params.Flag1 = 'ReadXYZ';

            params.OnsiteC1 = -0.39969;
            params.OnsiteC2 = -0.42011;
            params.OnsiteB = 2.3222;
            params.OnsiteN = -1.7795;

            params.cut = 24.0;

        case 'Slater_Koster_GBN'
            [value,Name] = textread('Input.in','%f%s');

            params.layer2 = value(2);
            params.layer1 = value(3);
            
            params.aG = value(4);
            params.aBN = value(5);
            
            params.theta = value(6);
            params.phi = value(7);

            params.eps = params.aG/params.aBN-1;
            params.eps2 = params.aBN/params.aG-1;

            params.CAB = 1.987/1000;
            params.aG = 3.5*pi/180;
            params.CAB_BN = 4.418/1000;
            params.phiAB_BN = 26.10*pi/180;
            params.phase_mode = 'Nofull';

            params.Flag1 = 'ReadXYZ';

            params.OnsiteC1 = -0.39969;
            params.OnsiteC2 = -0.42011;
            params.OnsiteB = 2.3222;
            params.OnsiteN = -1.7795;

            params.cut = 24.0;
    end
end
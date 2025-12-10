%function [Hamiltonian, neighbor_cell, disp_vec]=BuildHamiltonianGBN(kx,ky,pos,posc,CellLattice,Hsize,cut,params)
function [Hamiltonian, neighbor_cell, disp_vec] = BuildHamiltonianGBN(a1, a2, atom_position, atom_position_Init, N, neighbor_cell, params)
%function Ham = Hamilton_GBN(kx, ky, pos, posc, Cell, Hsize, cut)
%global OnsiteC1 OnsiteC2 OnsiteB OnsiteN aG aBN theta phi eps eps2 CAB phiAB CAB_BN phiAB_BN phase_mode

aG = params.aG;
aBN = params.aBN;
theta = params.theta;
phi = params.phi;
eps = params.eps;
eps2 = params.eps2;
CAB = params.CAB;
phiAB = params.aG;
CAB_BN = params.CAB_BN;
phiAB_BN = params.phiAB_BN;
phase_mode = params.phase_mode;
cut = params.cut;

TempSize1 = size(atom_position);
TempSize2 = size(atom_position_Init);
TempSize1 = TempSize1(1);
TempSize2 = TempSize2(1);

if (TempSize2 ~= TempSize1) || (TempSize2 ~= N)
    error('Size of Hamiltonian have mistakes')
end

Hsize = N;

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
x_v = atom_position(:, 1);  % x坐标
y_v = atom_position(:, 2);  % y坐标
z_v = atom_position(:, 3);  % z坐标
layer = atom_position(:, 4);  % 层数
sublattice = atom_position(:, 5);  % 子晶格
atom_type = atom_position(:, 6);  % 原子种类
dx_v = atom_position(:, 7);  % 原子种类
dy_v = atom_position(:, 8);  % 原子种类
d_v = atom_position(:, 9);  % 原子种类


% ================== 构建哈密顿量 ==================
% 初始化哈密顿量张量
%Ham = zeros(Hsize);
HamOnsite = zeros(Hsize);
HamHopping = zeros(Hsize);
Hamiltonian = zeros(Hsize,Hsize,N_neighbor);
%HAB = zeros(Hsize);

temp_disx = zeros(N_neighbor,1);    % These three values are for temp calculating nearest nine neighbor cell
temp_disy = zeros(N_neighbor,1);    
temp_dis = zeros(N_neighbor,1);

alpha = eps+1;
alpha2 = eps2+1;
theta_rad = theta*pi/180;
phirad = (phi)*pi/180;
rotMat = [cos(phirad),sin(phirad);-sin(phirad),cos(phirad)];
%% Diagonal Elements
for i = 1:Hsize
    if sublattice(i) == 3
        x = dx_v(i)*rotMat(1,1)+dy_v(i)*rotMat(1,2);
        y = dx_v(i)*rotMat(2,1)+dy_v(i)*rotMat(2,2);
        dx = x;
        dx = dx/aG;
        dy = y;
        dy = dy/aG;  
    else
        x = dx_v(i)*rotMat(1,1)+dy_v(i)*rotMat(1,2);
        y = dx_v(i)*rotMat(2,1)+dy_v(i)*rotMat(2,2);
        dx = -x;
        dx = dx/aBN;
        dy = -y;
        dy = dy/aBN;  
    end
    if sublattice(i) == 3
        if atom_type(i) == 1
            B = -0.39969;
            C = -0.4277;
            A = -0.37776;
            HamOnsite(i,i) = HamOnsite(i,i)+Interpolation(dx,dy,A,B,C,1);
        elseif atom_type(i) == 2
            B = -0.42011;
            C = -0.42048;
            A = -0.36739;
            HamOnsite(i,i) = HamOnsite(i,i)+Interpolation(dx,dy,A,B,C,1);
        else
            fprintf('Onsite, Carbon, sublattice have mistake\n')
            pos(i,:)
            pause
        end
    elseif sublattice(i) == 1
        B = 2.3222;
        C = 2.3476;
        A = 2.2731;
        HamOnsite(i,i) = HamOnsite(i,i)+Interpolation(dx,dy,A,B,C,1);
        if layer(i) == 1
            fprintf('B or N atom on bottom layer, check configuration\n')
        end
    elseif sublattice(i) == 2
        B = -1.7795;
        C = -1.7778;
        A = -1.8324;
        HamOnsite(i,i) = HamOnsite(i,i)+Interpolation(dx,dy,A,B,C,1);
        if layer(i) == 1
            fprintf('B or N atom on bottom layer, check configuration\n')
        end    
    end
end

Hamiltonian(:,:,1) = HamOnsite;

%% Hopping Elements
for i = 1:Hsize
    if sublattice(i) == 3
        x = dx_v(i)*rotMat(1,1)+dy_v(i)*rotMat(1,2);
        y = dx_v(i)*rotMat(2,1)+dy_v(i)*rotMat(2,2);
        dx = x;
        dx = dx/aG;
        dy = y;
        dy = dy/aG;  
    else
        x = dx_v(i)*rotMat(1,1)+dy_v(i)*rotMat(1,2);
        y = dx_v(i)*rotMat(2,1)+dy_v(i)*rotMat(2,2);
        dx = -x;
        dx = dx/aBN;
        dy = -y;
        dy = dy/aBN;  
    end
    G1 = 4.0*pi/(sqrt(3.0));
    HAB_ii = 2*CAB*cos(sqrt(3)*G1*dx/2)*cos(G1*dy/2-phiAB)-2*CAB*cos(G1*dy+phiAB)-1j*2*sqrt(3)*CAB*sin(sqrt(3)*G1*dx/2)*sin(G1*dy/2-phiAB);
    A_ii = real(HAB_ii);
    B_ii = imag(HAB_ii);
    delta1 = 2*A_ii/3;
    delta2 = -(A_ii+sqrt(3)*B_ii)/3;
    delta3 = (-A_ii+sqrt(3)*B_ii)/3;
    HAB_ii = 2*CAB_BN*cos(sqrt(3)*G1*dx/2)*cos(G1*dy/2-phiAB_BN)-2*CAB_BN*cos(G1*dy+phiAB_BN)-1j*2*sqrt(3)*CAB_BN*sin(sqrt(3)*G1*dx/2)*sin(G1*dy/2-phiAB_BN);
    A_ii_BN = real(HAB_ii);
    B_ii_BN = imag(HAB_ii);
    delta1_BN = 2*A_ii_BN/3;
    delta2_BN = -(A_ii_BN+sqrt(3)*B_ii_BN)/3;
    delta3_BN = (-A_ii_BN+sqrt(3)*B_ii_BN)/3;
    for j = 1:Hsize
        %[temp_disx,temp_disy,temp_dis] = periodic_condition(Cell,temp_disx,temp_disy,temp_dis,pos(i,:),pos(j,:),a1,a2);
        %temp_disz = z_v(i)-z_v(i);
        for k = 1:N_neighbor
            dx0 = disp_vec(k, 1);   % 当前近邻的x位移
            dy0 = disp_vec(k, 2);   % 当前近邻的y位移

            temp_disx(k) = x_v(i)-x_v(j)+dx0;        
            temp_disy(k) = y_v(i)-y_v(j)+dy0;
            temp_disz = z_v(i)-z_v(j);
            temp_dis(k) = sqrt(temp_disx(k)^2+temp_disy(k)^2);

            if temp_dis(k) < cut && ((temp_disx(k) ~= 0) || (temp_disy(k) ~= 0) || (temp_disz ~= 0))
                if strcmp(phase_mode, 'full')
                    exp_factor = exp(1i*(kx*temp_disx(k)+ky*temp_disy(k)));
                else
                    exp_factor = 1; % 幅度模式
                end
                if layer(i) == layer(j) && sublattice(i) == 3 && sublattice(j) == 3                    
                    R = [temp_disx(k),temp_disy(k),temp_disz,temp_dis(k)];
                    r0 = 2.438977765130281/sqrt(3);
                    HamHopping(i,j) = HamHopping(i,j)+Full_C(R,atom_position(i,:),atom_position(j,:),aG)*exp_factor;
                    Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + F2G2_C(R,atom_position(i,:),atom_position(j,:),aG)*exp_factor;
                    a_cc = aG/sqrt(3);
                    if atom_type(i) == 1 && abs(temp_disx(k)) < aG/2-0.001 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta1*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta1*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 1 && temp_disx(k) > 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta2*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta2*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 1 && temp_disx(k) < 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta3*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta3*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 2 && abs(temp_disx(k)) < aG/2-0.001 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta1*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta1*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 2 && temp_disx(k) > 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta3*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta3*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 2 && temp_disx(k) < 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta2*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta2*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        fprintf('Might have mistake\n')
                    end
                elseif layer(i) == layer(j) && not(sublattice(i) == 3 && sublattice(j) == 3)
                    R = [temp_disx(k),temp_disy(k),temp_disz,temp_dis(k)];
                    r0 = 2.438977765130281/sqrt(3);
                    HamHopping(i,j) = HamHopping(i,j)+F2G2_BN(R,atom_position(i,:),atom_position(j,:),aBN)*exp_factor;
                    Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + Full_BN(R,atom_position(i,:),atom_position(j,:),aBN)*exp_factor;
                    a_cc = aBN/sqrt(3);
                    if atom_type(i) == 3 && abs(temp_disx(k)) < aG/2-0.001 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta1_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta1_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 3 && temp_disx(k) > 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta2_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta2_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 3 && temp_disx(k) < 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta3_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta3_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 4 && abs(temp_disx(k)) < aG/2-0.001 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta1_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta1_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 4 && temp_disx(k) > 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta3_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta3_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif atom_type(i) == 4 && temp_disx(k) < 0 && temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        HamHopping(i,j) = HamHopping(i,j)+delta2_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                        Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + delta2_BN*exp_factor*exp(-3.37*(R(4)-r0)/r0);
                    elseif temp_dis(k) < a_cc+0.1*a_cc && temp_dis(k) > a_cc-0.1*a_cc
                        fprintf('Might have mistake\n')
                    end
                elseif layer(i) ~= layer(j)
                    R = [temp_disx(k),temp_disy(k),temp_disz,temp_dis(k)];
                    HamHopping(i,j) = HamHopping(i,j)+Interlayer(R,atom_position(i,:),atom_position(j,:))*exp_factor;
                    Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + Interlayer(R,atom_position(i,:),atom_position(j,:))*exp_factor;
                else
                    fprintf('Some near neighbor is not defined in the Hamiltonian')
                    pos(i,:)
                    pos(j,:)
                    HamHopping(i,j) = HamHopping(i,j)+0;
                    Hamiltonian(i, j, k) = Hamiltonian(i, j, k) + 0;
                end
            end
        end
    end
end

%Ham =  (HamHopping+HamHopping')/2+HamOnsite;



    

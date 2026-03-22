%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TII文章-面向可靠性的传感器调度
% 不同的故障分布的对比图
% 不仅要画出来求解的时间，也要画出来优化的效果
%
% Peizhe Li
% 2026-01-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Paramater setting
% --------------------------
% 参数定义（物理对象参数）
% --------------------------
c = 460;             % 比热容 J/(kg·K)
rho = 7.85e3;        % 密度 kg/m^3
lambda = 40;         % 热导率 W/(m·K)
delta_b = 0.05;      % 每个厚度单元的厚度 (m)
delta_l = 5;         % 每个长度单元长度 (m)
delta_t = 0.1;       % 时间步长 (s)
v = 0/60;               % 输送速度 (m/s)
% v=0
epsilon = 0.85;
sigma0 = 5.67e-8;
x_inf = 60 + 273.15;     % 恒定外界温度
x0 = 1180;        % 初始温度为1180摄氏度
T_total = 200/delta_t;   % 仿真总步数

% --------------------------
% 钢板空间离散结构
% --------------------------
nu = 3;                % 厚度方向分段数
tau_s = 5;            % 长度方向分段数
n = nu * tau_s;        % 状态变量维度

% --------------------------
% 传感器数量设定并初始化运行时间参数
% --------------------------
Sensor_total = 10;
ite = 1;

% --------------------------
% 优化算法运行参数
% --------------------------
T_Simulation = 2000;
gamma_hist = cell(length(Sensor_total),1);
L_mat_hist = cell(length(Sensor_total),1);
time_hist = zeros(1,length(Sensor_total));
Upsilon = 10;

num_sensors = Sensor_total;
num_subproblem = 20;
time_sub_avg = zeros(length(Sensor_total),10);
time_opt_avg = zeros(1,10);
max_ite = 1;
time_sub_sum = [];
time_opt_sum = [];

W = 0.1 * eye(n);       % 过程噪声协方差
R = 0.1 * eye(num_sensors);  % 观测噪声协方差

% --------------------------
% 优化问题参数
% --------------------------
zeta = 15;
% --------------------------
% 通信信道参数
% --------------------------
% N0 = 10^((-75-30)/10);  % 信道噪声
% Bit = 64*8;             % 数据包大小
% B = 2e6;                % 带宽
% [distances, gains] = compute_channel_gain(num_sensors);
% transmit_power = (N0 * B) ./ gains .* (2.^(Bit / B / 0.01) - 1);
%

% 参数设置
slot = 5;
N0 = 10^((-75-30)/10);  % 信道噪声
transmit_power = 10^((15-30)/10);
[distances, gains] = compute_channel_gain(num_sensors);
SNR0 = 500;     % SNR门限
% p_i = 10;     % 发射功率
% 理论概率
x = (N0 * SNR0 ) ./ (transmit_power * gains);
P_success = exp(-x); % 每个传感器的传输成功率

% --------------------------
% 设备温度模型参数
% 传感器可靠性参数，看xiaofan yu开源的代码确认一下
% --------------------------
aT = 0.9;
bT = 0.1;
cT = 0.03;
Tm = 60 + 273.15;

pu = 15*10^(-3);

%% 构建系统模型
% --------------------------
% 差分矩阵 Q（厚度方向的二阶中心差分）
% --------------------------
Q = zeros(nu);
Q(1,1) = -2; Q(1,2) = 2;
for i = 2:nu-1
    Q(i,i-1) = 1;
    Q(i,i)   = -2;
    Q(i,i+1) = 1;
end
Q(nu,nu-1) = 2;
Q(nu,nu) = -2;

% --------------------------
% Lambda 矩阵（长度方向的对流项）
% --------------------------
% omega = v / (2 * delta_l);
omega = delta_t*v / (delta_l);
Lambda = zeros(n);
I_nu = eye(nu);

for i = 1:tau_s
    idx = (i-1)*nu + 1 : i*nu;
    Lambda(idx, idx) = Lambda(idx, idx) + (1-omega) * I_nu;
    if i < tau_s
        Lambda(idx, idx+nu) = Lambda(idx, idx+nu) + omega * I_nu;
    end
end

% --------------------------
% D_j 矩阵（热扩散系数）
% --------------------------
% alpha = delta_t * lambda / (delta_b^2 * rho * c);
D_j = delta_t * lambda / (delta_b^2 * rho * c) * eye(nu);    % 局部块的热扩散项

% --------------------------
% 状态转移矩阵 F 构造
% --------------------------
F = Lambda';
for i = 1:tau_s
    idx = (i-1)*nu + 1 : i*nu;
    F(idx, idx) = F(idx, idx) + D_j * Q;
end

% --------------------------
% 构建 B 矩阵
% --------------------------
B_coeff = 2 * epsilon * sigma0 * delta_t / (rho * c * delta_b);
for i = 1:tau_s
    idx_top = (i - 1) * nu + 1;        % 每段顶部索引
    idx_bottom = i * nu;              % 每段底部索引
    B(idx_top,1) = B_coeff;
    B(idx_bottom,1) = B_coeff;
end

% --------------------------
% 构建 C 矩阵
% --------------------------
sensor_positions = round(linspace(1, tau_s, num_sensors));  % 观测的 slab 编号（长度方向）
% 随机指定每个传感器观测上表面（1）或下表面（nu）
surface_choice = randi([1, 2], num_sensors, 1);  % 1 表示上表面，2 表示下表面
C = zeros(num_sensors, n);
for i = 1:num_sensors
    slab_idx = sensor_positions(i);
    if surface_choice(i) == 1
        state_idx = (slab_idx - 1) * nu + 1;  % 上表面
    else
        state_idx = slab_idx * nu;           % 下表面
    end
    C(i, state_idx) = 1;
end

%% 优化算法流程，从不同数量的子问题和不同规模的传感器数量去分析
% 和其他调度策略对比可靠性，事件触发，最小化能耗，最小化AoI；对比估计的精度和可靠性，先写出来完整的优化问题
% observable_sensors = find_observable_sensors_per_jordan_block(F, C);    % 把每个jordan块对应能观测的传感器找到
[Gth, eigvals] = extract_Gth_jordan(F, C);  % 计算每个jordan块对应的传感器

L_eigs = run_algo1_per_eigenvalue(F, C, Gth);
L_eigs_all = merge_and_deduplicate_L_eigs(L_eigs);

L_all = cell(1);
L_all{1} = L_eigs{1};
L_mat = L_all{1};  % 这里得到了整体的L矩阵

N_M = 0;
for i = length(L_all)
    N_M = N_M + length(L_all{i}); % 这里是得到了整体的N_L是多大
end

% 初始化温度就等于环境温度
T_init = Tm*ones(num_sensors,1);
p = 0.05 + transmit_power;

%% 对subproblem数量进行循环
gamma1 = [];T_init_rec=[];
gamma1 = [];T_init_rec=[];
T_init = Tm*ones(num_sensors,1);
time_mean = [];
for tq = 1:Upsilon:T_Simulation
    num_subproblem = factor(N_M);
    num_subproblem = num_subproblem(end);

    X         = cell(num_subproblem,1);
    Obj_value = cell(num_subproblem,1);
    A_total   = cell(num_subproblem,1);
    b_total   = cell(num_subproblem,1);
    L_local   = cell(num_subproblem,1);

    A_sub = []; b_sub = []; L_sub = [];

    start_point_A = zeros(num_subproblem,1);
    end_point_A   = zeros(num_subproblem,1);
    start_point_L = zeros(num_subproblem,1);
    end_point_L   = zeros(num_subproblem,1);

    % N_mod = mod(N_M,num_subproblem);
    % L_mat = [L_mat;L_mat(1:num_subproblem-N_mod,:)];
    N_M = size(L_mat,1); N_mod = mod(N_M,num_subproblem);
    end_point = 0;

    for i = 1:num_subproblem
        L_local{i} = L_mat((i-1)*(N_M/num_subproblem)+1:i*(N_M/num_subproblem),:);
        [A_total{i},b_total{i}] = local_constraint(L_local{i}, num_sensors, Upsilon, slot, P_success); % 通过这个函数实现对约束的定义
        start_point_A(i,1) = size(A_sub,1)+1;
        start_point_L(i,1) = size(L_sub,1)+1;
        A_sub = [A_sub; A_total{i}];
        b_sub = [b_sub; b_total{i}];
        L_sub = [L_sub; L_local{i}];
        end_point_A(i,1) = size(A_sub,1);
        end_point_L(i,1) = size(L_sub,1);
    end

    % 开最大数量的并行
    time_sub = zeros(1,num_subproblem);
    parfor step = 1:num_subproblem
        Local_data = L_sub(start_point_L(step,1):end_point_L(step,1),:);
        Local_A = A_sub(start_point_A(step,1):end_point_A(step,1),:);
        Local_b = b_sub(start_point_A(step,1):end_point_A(step,1),:);
        [X_value, f_value, time] = slove_X(Local_data, p, zeta, T_init, aT, bT, cT, Tm, num_sensors, Upsilon, Local_A, Local_b);
        X{step}= X_value; Obj_value{step} = f_value; time_sub(1,step) = time;
    end
    time_mean = [time_mean time_sub];
    % 输出优化结果
    % 首先把 Obj_value 转换为数值数组
    Obj_numeric = cellfun(@(c) c, Obj_value);

    % 找到最小值的位置
    [~, idx_min] = min(Obj_numeric);

    % 输出对应的 X
    X_min = X{idx_min};
    gamma1 = [gamma1 reshape(X{idx_min}(1:num_sensors*Upsilon,end), Upsilon, num_sensors)'];

    % 更新对应的初始温度
    gamma_now = reshape(X{idx_min}(1:num_sensors*Upsilon,end), Upsilon, num_sensors)';
    for i = 1:num_sensors
        T_init(i,1) = compute_temperature_after_upsilon(gamma_now(i,:), T_init(i,1), p, aT, bT, cT, Tm);
    end
    T_init_rec=[T_init_rec T_init];
end
gamma_hist{ite} = gamma1;
time_hist(1,ite) = mean(time_mean);
L_mat_hist{ite} = L_mat;
ite = ite+1;


%% 画出相应的对比图
% 首先应该对比估计性能
% gamma1 = gamma;
p = 0.05 + transmit_power;
power = p;
Gain_Simulation = 20000/T_Simulation;
gamma = cell(length(Sensor_total),1); Tem = cell(length(Sensor_total),1);
C1 = cell(length(Sensor_total),1); P = cell(length(Sensor_total),1);
K = cell(length(Sensor_total),1); P_trace = cell(length(Sensor_total),1);

for i = 1:length(Sensor_total)
    gamma{i} = kron(ones(1,Gain_Simulation),gamma_hist{i});
    P{i} = 30*eye(size(F,1));
    P_trace{i} = zeros(1,T_Simulation*Gain_Simulation);
    Tem{i} = zeros(num_sensors,T_Simulation*Gain_Simulation);
    Tem{i}(:,1) = Tm*ones(num_sensors,1);
end


for i = 1:T_Simulation*Gain_Simulation
    % 对每个结果都做一遍这个循环
    for j = 1:length(Sensor_total)
        % 把对应的观测矩阵拿出来
        idx = find(gamma{j}(:,i) ~= 0);
        C1{j} = C(idx,:);

        % 协方差的一步预测
        P{j} = F*P{j}*F'+W;

        % 协方差的更新
        K{j} = P{j}*C1{j}'/(C1{j}*P{j}*C1{j}'+0.1*eye(size(C1{j},1)));
        P{j} = (eye(size(K{j},1))-K{j}*C1{j})*P{j};

        % 协方差的迹存下来
        P_trace{j}(1,i) = trace(P{j});
        % 每个传感器温度的更新
        for k = 1:Sensor_total(j)
            Tem{j}(k,i+1) = compute_temperature_after_upsilon(gamma{j}(k,i), Tem{j}(k,i), power, aT, bT, cT, Tm);
        end
    end
end

% 直接用平均温度去算可靠性增益
Ea_total = 0.7:0.05:1.0;
Tem_avg = cell(length(Ea_total),1);
MTTF = cell(length(Ea_total),1);

for j = 1:length(Ea_total)
    Tem_avg{j} = mean(Tem{1},2);
end
sensor_params = repmat(struct(), 1, num_sensors);
rng(130);
for j = 1:length(Ea_total)
    for i = 1:num_sensors
        sensor_params(i).delta = 2;
        sensor_params(i).kappa = 1.5e1 + 1*randn();
        sensor_params(i).J = 1.5e6 + 100*randn();
        % sensor_params(i).kappa = 6.5e2;
        % sensor_params(i).J = 1.5e6;
        sensor_params(i).Jc = 2.5e5;
        sensor_params(i).Ea = Ea_total(j);
        sensor_params(i).B0 = 8.617e-5;
        sensor_params(i).aT = 0.9;
        sensor_params(i).bT = 0.15;
        sensor_params(i).cT = 0.03;
        sensor_params(i).Tm = 60 + 273.15;
        sensor_params(i).pi_t = 0.015 + transmit_power;
        p = sensor_params(i);
        MTTF{j}(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem_avg{j}(i))/24/30/12;
    end
end

[N_L,g_l] = size(L_mat_hist{1});
for k = 1:length(Ea_total)
    MTTF_event = 10^2000*ones(length(Ea_total),N_L);
    % 直接利用L_mat去计算系统的MTTF
    for i = 1:N_L
        for j = 1:g_l
            MTTF_event(k,i) = min([MTTF_event(k,i),MTTF{k}(L_mat_hist{1}(i,j))]);
        end
    end
    MTTF_all(k) = max(MTTF_event(k,:));

    % 还需要对比一下整体的能耗
    Energy(k) = sum(sum(gamma{1}))*transmit_power*delta_t;
end

% plot result
% 先对比估计效果
color_plot = [
    0.929, 0.694, 0.125;  % 黄色
    0.200, 0.600, 0.800;  % 蓝
    0.200, 0.800, 0.400;  % 绿
    0.850, 0.500, 0.000;  % 棕
    0.800, 0.300, 0.400;  % 玫红
    0.900, 0.600, 0.300;  % 珊瑚
    ];
figure(2)
for i = 1:length(Sensor_total)
    hold on;
    plot([1:T_Simulation*Gain_Simulation]*0.1, P_trace{i}, '-.', 'Color', color_plot(i,:), ...
        'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', ['N=' num2str(Sensor_total(i))]);
end
% 坐标和标签
xlabel('Time');
ylabel('tr(P)');
legend('Location', 'northwest');
% title('Comparison of Solving Time')
grid on;
box on;
set(gca, 'FontSize', 15);
% ylim([0,20])

% 画图，整理成数组的形式，由于选择的weight不同，把对应的总cost的不要了
Bar_Data = [];
for i = 1:length(Ea_total)
    Bar_Data = [Bar_Data;MTTF_all(i)/MTTF_all(1)];
end
% Bar_Lab  = categorical({'$\Upsilon=4$','$\Upsilon=6$','$\Upsilon=8$','$\Upsilon=10$'});
% Bar_Lab  = reordercats(Bar_Lab,{'$\Upsilon=4$','$\Upsilon=6$','$\Upsilon=8$','$\Upsilon=10$'});
% 自动生成标签
LabelStr = arrayfun(@(x) sprintf('$%.2f$', x), Ea_total, 'UniformOutput', false);
% 创建和排序categorical
Bar_Lab = categorical(LabelStr);
Bar_Lab = reordercats(Bar_Lab, LabelStr);
figure(6)
Cost_bar = bar(Bar_Lab,Bar_Data);
Cost_bar(1).FaceColor = '#7fb8dd';
% Cost_bar(2).FaceColor = '#eba88b';
% Cost_bar(3).FaceColor = '#f5d78f';
ylim([0.5, 10^20]);
xlabel('Activation energy (eV)','Interpreter','latex')
% 显式替换 x 轴标签为 LaTeX 公式
xticklabels(LabelStr)
set(gca, 'TickLabelInterpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',16)
ylabel('Relative MTTF')
set(gca, 'YScale', 'log');

%% 不同的电流密度下的可靠性
sensor_params = repmat(struct(), 1, num_sensors);
J_total = 1.0:0.1:1.6;
rng(130);
for j = 1:length(J_total)
    for i = 1:num_sensors
        sensor_params(i).delta = 2;
        sensor_params(i).kappa = 1.5e1 + 1*randn();
        sensor_params(i).J = J_total(j)*10^6 + 100*randn();
        % sensor_params(i).kappa = 6.5e2;
        % sensor_params(i).J = 1.5e6;
        sensor_params(i).Jc = 2.5e5;
        sensor_params(i).Ea = 0.9;
        sensor_params(i).B0 = 8.617e-5;
        sensor_params(i).aT = 0.9;
        sensor_params(i).bT = 0.15;
        sensor_params(i).cT = 0.03;
        sensor_params(i).Tm = 60 + 273.15;
        sensor_params(i).pi_t = 0.015 + transmit_power;
        p = sensor_params(i);
        MTTF{j}(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem_avg{j}(i))/24/30/12;
    end
end

[N_L,g_l] = size(L_mat_hist{1});
for k = 1:length(J_total)
    MTTF_event = 10^2000*ones(length(Ea_total),N_L);
    % 直接利用L_mat去计算系统的MTTF
    for i = 1:N_L
        for j = 1:g_l
            MTTF_event(k,i) = min([MTTF_event(k,i),MTTF{k}(L_mat_hist{1}(i,j))]);
        end
    end
    MTTF_all(k) = max(MTTF_event(k,:));

    % 还需要对比一下整体的能耗
    Energy(k) = sum(sum(gamma{1}))*transmit_power*delta_t;
end

% plot result
% 画图，整理成数组的形式，由于选择的weight不同，把对应的总cost的不要了
Bar_Data = [];
for i = 1:length(J_total)
    Bar_Data = [Bar_Data;MTTF_all(i)/MTTF_all(1)];
end
% Bar_Lab  = categorical({'$\Upsilon=4$','$\Upsilon=6$','$\Upsilon=8$','$\Upsilon=10$'});
% Bar_Lab  = reordercats(Bar_Lab,{'$\Upsilon=4$','$\Upsilon=6$','$\Upsilon=8$','$\Upsilon=10$'});
% 自动生成标签
LabelStr = arrayfun(@(x) sprintf('$%.1f$', x), J_total, 'UniformOutput', false);
% 创建和排序categorical
Bar_Lab = categorical(LabelStr);
Bar_Lab = reordercats(Bar_Lab, LabelStr);
figure(6)
Cost_bar = bar(Bar_Lab,Bar_Data);
Cost_bar(1).FaceColor = '#7fb8dd';
% Cost_bar(2).FaceColor = '#eba88b';
% Cost_bar(3).FaceColor = '#f5d78f';
ylim([0, 1.2]);
xlabel('Current density ($\times 10^6$ A/cm$^2$)','Interpreter','latex')
% 显式替换 x 轴标签为 LaTeX 公式
xticklabels(LabelStr)
set(gca, 'TickLabelInterpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',16)
ylabel('Relative MTTF')
% set(gca, 'YScale', 'log');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 函数定义
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 能观性检查
function Noise = Gaussian_Noise_Generate(Cov,T)
R     = chol(Cov);
[m,~] = size(R);
Noise = R*randn(m,T);
end

function [isObservable, ObsvMatrix, rankObsv, nStates] = checkObservability(A, C)
% CHECKOBSERVABILITY 检查线性系统的能观性
%   输入:
%       A - 系统矩阵 (n x n)
%       C - 观测矩阵 (p x n)
%   输出:
%       isObservable - 布尔值，true表示系统能观
%       ObsvMatrix   - 能观性矩阵
%       rankObsv     - 能观性矩阵的秩
%       nStates      - 系统状态数

% 验证输入矩阵维度
[n, m] = size(A);
[p, q] = size(C);
if n ~= m
    error('A必须为方阵 (n x n)');
end
if q ~= n
    error('C的列数必须与A的行数一致');
end

nStates = n;  % 系统状态数
ObsvMatrix = []; % 初始化能观性矩阵

% 手动构造能观性矩阵: [C; C*A; C*A^2; ... ; C*A^(n-1)]
for k = 0:nStates-1
    ObsvMatrix = [ObsvMatrix; C * (A^k)];
end

% 计算能观性矩阵的秩
rankObsv = rank(ObsvMatrix);

% 判断是否满秩 (能观性条件)
isObservable = (rankObsv >= nStates); % 使用 >= 以应对数值计算误差
end

% 计算钢板参数
function [rho, lambda_, c, epsilon, sigma] = steel_properties(x)
% 输入: x = 温度 (单位 K)
% 输出: rho, lambda_, c, epsilon, sigma0（热辐射常数）

% --- 密度 rho [kg/m^3] ---
rho = 7843.76 - 0.2958 * (x - 273) - 5.65e-5 * (x - 273).^2;

% --- 导热率 lambda [W/mK] ---
lambda_ = 20.14 + 9.313e-3 * (x - 273);

% --- 比热容 c [J/kgK]（分段函数） ---
c = zeros(size(x));
for i = 1:length(x)
    xi = x(i)-273;
    if xi >= 800 && xi < 1000
        c(i) = 4.583*xi - 4720.3 + 1.109e9 / xi^2;
    elseif xi >= 1000 && xi < 1042
        c(i) = 12.476 * xi - 11501;
    elseif xi >= 1042 && xi < 1060
        c(i) = -32 * xi + 34871.2;
    elseif xi >= 1060 && xi < 1184
        c(i) = 5.987 * xi - 10068.18 + 5.21e9 / xi^2;
    elseif xi >= 1184 && xi <= 1665
        c(i) = 0.15 * xi + 429.85;
    else
        c(i) = NaN; % 超出定义范围
    end
end

% --- 发射率 epsilon(x) ---
epsilon = ((x - 273) / 1000) .* (0.125 * (x - 273) / 1000 - 0.38) + 1.1;

% --- 热辐射常数 sigma [W/m^2K^4] ---
sigma = 5.67e-8;  % 常量

end

% 传感器到远程估计器的距离和信道增益
function [distances, gains] = compute_channel_gain(N)
% 输入：N = 传感器数量
% 输出：distances = 每个传感器到估计器的距离
%         gains = 每个传感器的信道增益

dx = 10;                  % 传感器间隔
lossfactor = 1;            % 信道损耗因子
sensor_x = (0:N-1) * dx;   % 横向坐标，假设第一个传感器在0
sensor_y = zeros(1, N);    % 所有传感器纵向坐标为0

% 远程估计器坐标
estimator_x = mean(sensor_x);  % 横向在中点
estimator_y = 7;               % 纵向固定7米

% 计算距离
distances = sqrt((sensor_x - estimator_x).^2 + (sensor_y - estimator_y).^2);

% 计算信道增益
gains = lossfactor * distances.^(-3);
end

% 计算平均温度
function T_avg = compute_avg_temperature(aT, bT, cT, Tm, pi_t, T0)
% 输入：温度模型参数与功耗序列
% 输出：平均温度（标量）

if nargin < 6
    T0 = Tm;  % 初始温度默认等于环境温度
end

n = 10^5;
T = zeros(1, n);
T(1) = T0;

for k = 2:n
    T(k) = aT * T(k-1) + bT * pi_t + cT * Tm;
end

T_avg = mean(T);
end


% 计算设备MTTF
function MTTF = estimate_mttf_from_Tavg(delta, kappa, J, Jc, Ea, B0, T_avg)
% 输入：Weibull参数与平均温度
% 输出：近似MTTF

gamma_term = gamma(1 + 1/delta);
mu = (kappa * (J - Jc)^(-2)) * exp(Ea / (B0 * T_avg)) / gamma_term;

MTTF = mu * gamma_term;  % 近似：MTTF ≈ μ * Gamma
end

% 计算每个jordan块对应的传感器
function [Gth, eigvals] = extract_Gth_jordan(F, C)
% 输出：
% Gth{u}{j} 表示第 u 个特征值的第 j 个 Jordan 块对应的可观测传感器
% eigvals(u) 是第 u 个特征值

[V, J] = jordan(F);               % Jordan 分解
Cj = C * V;                       % 变换到 Jordan 坐标下
n = size(F, 1);
lambda_list = diag(J);

% 提取所有不重复的特征值
eigvals = unique(lambda_list, 'stable');
Gth = cell(1, length(eigvals));

% 遍历 Jordan 块
i = 1;
while i <= n
    lambda = J(i, i);
    % 找该特征值在 eigvals 中的索引
    u = find(abs(eigvals - lambda) < 1e-8, 1);

    % 当前块的维度（行内连续λ，超对角为1）
    blk_size = 1;
    while i + blk_size <= n && J(i+blk_size, i+blk_size) == lambda ...
            && abs(J(i+blk_size-1, i+blk_size) - 1) < 1e-8
        blk_size = blk_size + 1;
    end

    blk_start = i;  % 该 Jordan 块的首状态索引
    observable_sensors = [];

    % 遍历每个传感器（Cj的行）
    for s = 1:size(Cj, 1)
        if abs(Cj(s, blk_start)) > 1e-8
            observable_sensors(end+1) = s;
        end
    end

    % 记录该 Jordan 块的观测传感器编号
    if isempty(Gth{u})
        Gth{u} = {};
    end
    Gth{u}{end+1} = observable_sensors;

    % 移动到下一个 Jordan 块
    i = i + blk_size;
end
end

% 返回矩阵的代数重数
function geom_mult_matrix = jordan_geometric_multiplicity(F)
[~, J] = jordan(F);
n = size(J,1);
lambda_diag = diag(J);

unique_lambdas = unique(lambda_diag, 'stable');  % 按出现顺序提取不重复特征值
geo_mults = zeros(size(unique_lambdas));

for i = 1:length(unique_lambdas)
    lambda = unique_lambdas(i);
    count = 0;

    for row = 1:n
        % 判断当前是 λ 的 Jordan 块起始：
        if J(row, row) == lambda
            if row == 1 || J(row-1, row) ~= 1
                count = count + 1;
            end
        end
    end
    geo_mults(i) = count;
end

geom_mult_matrix = [unique_lambdas.'; geo_mults.'];
end

% 构造\bar{C}
function [C_bar, blk_info] = construct_C_bar(C, A)
[V, J] = jordan(A);
Cj = C * V;
n = size(A,1);

blk_starts = [];
blk_lambdas = [];

i = 1;
while i <= n
    lambda = J(i,i);
    blk_starts(end+1) = i;
    blk_lambdas(end+1) = lambda;

    blk_size = 1;
    while i+blk_size <= n && J(i+blk_size, i+blk_size) == lambda ...
            && abs(J(i+blk_size-1, i+blk_size)-1) < 1e-8
        blk_size = blk_size + 1;
    end
    i = i + blk_size;
end

% 排序 Jordan 块：按特征值分组
[sorted_lambdas, sort_idx] = sort(blk_lambdas);  % ascending order
sorted_blk_starts = blk_starts(sort_idx);

% 构造 C_bar
C_bar = zeros(size(C,1), length(sorted_blk_starts));
for k = 1:length(sorted_blk_starts)
    C_bar(:,k) = Cj(:, sorted_blk_starts(k));
end

% 构造 blk_info
blk_info.block_start = sorted_blk_starts;
blk_info.eigenvalue = sorted_lambdas;

% 构造 lambda_group_end: 每个特征值在C_bar中最后一个列索引
unique_lambdas = unique(sorted_lambdas, 'stable');
lambda_group_end = zeros(1, length(unique_lambdas));

cursor = 0;
for i = 1:length(unique_lambdas)
    count = sum(abs(sorted_lambdas - unique_lambdas(i)) < 1e-8);
    cursor = cursor + count;
    lambda_group_end(i) = cursor;
end

blk_info.lambda_group_end = lambda_group_end;
end

% 循环运行算法1，对每个特征值
function L_all_eigs = run_algo1_per_eigenvalue(F, C, Gth)
% 计算所有特征值的几何重数
geom_mult_mat = jordan_geometric_multiplicity(F);  % 2×k 矩阵
lambdas = geom_mult_mat(1, :);
geom_mults = geom_mult_mat(2, :);

% 存储每个 λ 的L结果：cell{j} 对应 λ_j
L_all_eigs = cell(1, length(lambdas));

[C_bar, blk_info] = construct_C_bar(C, F);

interval = [0,blk_info.lambda_group_end];
for j = 1:length(lambdas)
    gi = geom_mults(j);  % 当前特征值的几何重数
    L_all_eigs{j} = hierarchical_independent_sets(C_bar(:,interval(j)+1:interval(j+1)), Gth{j}, gi);

    fprintf('Eigenvalue λ = %.3f, geometric multiplicity gi = %d\n', ...
        lambdas(j), gi);
end
end

% 实现 Algorithm 1
function L_all = hierarchical_independent_sets(C, Gthu, gi)
% 输入：
%   - C:     观测矩阵（m × n）
%   - Gthu:  每一层的候选传感器集合，cell{1:gi}
%   - gi:    总层数（通常是特征值的几何重数）
% 输出：
%   - L_all: cell，L_all{i} 是第 i 个初始传感器对应的所有最大无关子集，行 × gi 矩阵

G1 = Gthu{1};
N = length(G1);
L_all = [];  % 每个初始传感器一个结果矩阵

for idx = 1:N
    i = G1(idx);  % 第 idx 个初始传感器编号

    % ---- 初始化 L^1_i ----
    L_prev = i;  % 第一层只有一个子集，就是{i}

    for tau = 2:gi
        L_curr = [];  % 当前 tau 层所有合法子集
        for k = 1:size(L_prev, 1)  % 遍历上一层的所有子集
            curr_subset = L_prev(k, :);
            C_prev = C(curr_subset, 1:tau);

            for j = Gthu{tau}
                if ismember(j, curr_subset)
                    continue;  % 避免重复
                end
                candidate = [curr_subset, j];
                C_try = C(candidate, 1:tau);
                if rank(C_try) > rank(C_prev)
                    L_curr = [L_curr; candidate];
                end
            end
        end
        if isempty(L_curr)
            break;  % 当前层构不出新子集，提前终止
        end
        L_prev = L_curr;  % 准备下一层
    end

    % 保存该初始传感器的所有最终路径（行 × gi 矩阵）
    L_all = [L_all;L_prev];
end
end

function L_eigs_all = merge_and_deduplicate_L_eigs(L_eigs)
% 拼接所有 cell 内的矩阵
L_combined = [];
for i = 1:length(L_eigs)
    if ~isempty(L_eigs{i})
        L_combined = [L_combined; L_eigs{i}];
    end
end

% 删除重复行
L_eigs_all = unique(L_combined, 'rows');
end

% 计算一段时间内传感器的平均温度
function T_final = compute_temperature_after_upsilon(gamma_i, T_init, p_i, a_T, b_T, c_T, T_m)
% gamma_i: 1 × Upsilon，传感器在 Upsilon 时间段内的调度激活序列（0/1）
% T_init: 初始温度 T_i(t_q)
% p_i:    传感器功耗
% a_T, b_T, c_T, T_m: 参数
%
% 返回：T_i(t_q + Upsilon)

T = T_init;
Upsilon = length(gamma_i);
for k = 1:Upsilon
    T = a_T * T + b_T * p_i * gamma_i(k) + c_T * T_m;
end
T_final = T;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 优化问题的目标函数和约束条件建模
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = fobj_P2(gamma, L_all, p, zeta, T_init, a_T, b_T, c_T, T_m)
% gamma: N × Υ, binary activation matrix
% p:     N × 1, sensor power weights
% T_init: N × 1, initial temperature T_i(t_q)
% a_T, b_T, c_T, T_m: scalars
[N, Upsilon] = size(gamma);
T_bar = zeros(N, 1);

for i = 1:N
    temp_sum = 0;
    for k = 0:Upsilon-1
        for j = 0:k-1
            temp_sum = temp_sum + (a_T^(k-1-j)) * (b_T*p(i) * gamma(i,j+1) + c_T + T_m);
        end
    end
    term1 = temp_sum / Upsilon;
    term2 = sum(a_T.^(0:Upsilon-1)) * T_init(i) / Upsilon;
    T_bar(i) = term1 + term2;
end

term1 = 0;
for l = 1:length(L_all)
    L_mat = L_all{l};       % size N_L × g_l
    N_L = size(L_mat, 1);
    g_l = size(L_mat, 2);

    sum_T = 0;
    for j = 1:N_L
        for i = 1:g_l
            idx = L_mat(j,i);         % sensor index
            sum_T = sum_T + T_bar(idx);
        end
    end
    term1 = term1 + sum_T / N_L;
end

term2 = sum(sum(p(:) .* gamma));  % energy cost
obj = term1 - zeta * term2;
end


%% initial value
function [X_value, f_value, time] = slove_X(L_mat , p, zeta, T_init, a_T, b_T, c_T, T_m, num_sensors, Upsilon, A_total, b_total)
% X: NΥ × 1, binary activation vector
% p:     N × 1, sensor power weights
% T_init: N × 1, initial temperature T_i(t_q)
% a_T, b_T, c_T, T_m: scalars

% 计算第1部分目标函数
N = num_sensors;

N_L = size(L_mat, 1);
g_l = size(L_mat, 2);

p = p*ones(N,1);

% X = binvar(num_sensor*Upsilon+size(L_mat , 1), 1);
X = binvar(size(A_total, 2)-1, 1);
b = sdpvar(1,1);

% 计算平均温度，同时也要计算出来0-1选择变量对应的约束
Aeq1 = []; beq1 = [];
Aeq2 = []; beq2 = [];
for i = 1:N
    temp_sum = 0;
    for k = 0:Upsilon-1
        for j = 0:k-1
            temp_sum = temp_sum + (a_T^(k-1-j)) * (b_T*p(i) * X((i-1)*Upsilon+1+j) + c_T + T_m);
        end
    end
    term1 = temp_sum / Upsilon;
    term2 = sum(a_T.^(0:Upsilon-1)) * T_init(i) / Upsilon;
    T_bar(i) = term1 + term2;
    Aeq1 = [Aeq1;zeros(1,num_sensors*Upsilon+(i-1)*(Upsilon+1)), ones(1,Upsilon+1),zeros(1,(N-i)*(Upsilon+1)+N_L)];
    beq1 = [beq1;ones(1)];
    Aeq2 = [Aeq2;zeros(1,(i-1)*Upsilon),ones(1,Upsilon),zeros(1,(N-i)*Upsilon), zeros(1,(i-1)*(Upsilon+1)), ...
        [0:-1:-Upsilon],zeros(1,(N-i)*(Upsilon+1)+N_L)];
    beq2 = [beq2;zeross(1)];
end

sum_T = 0;
for j = 1:N_L
    for i = 1:g_l
        idx = L_mat(j,i);         % sensor index
        sum_T = sum_T + T_bar(idx);
    end
end

term1 = sum_T / N_L;

term2 = sum(sum(kron(p',ones(1,Upsilon))* X(1:num_sensors*Upsilon,:)));  % energy cost
Objective = (term1 + zeta * term2);

Aeq3 = [zeros(1,num_sensors*Upsilon+num_sensors*(Upsilon+1)),ones(1,N_L)]; beq3 = 1;
Constraints = [A_total*[X;b] <= b_total;...
    [Aeq1;Aeq2;Aeq3]*X == [beq1;beq2;beq3];
    ];

options = sdpsettings('solver', 'gurobi', 'verbose', 1);
result = optimize(Constraints, Objective, options);

X_value = [value(X);value(b)];
f_value = value(Objective);
time = result.solvertime;

end

%% generation of constraint
function [A_total, b_total] = local_constraint(L_mat, num_sensors, Upsilon, slot, P_success)
% gamma 的排列顺序是每个按传感器来的，先是第一个传感器的Upsilon个调度决策，然后是第二个……

N_L = size(L_mat, 1);
g_l = size(L_mat, 2);

M = 10^3;

n = num_sensors*Upsilon+num_sensors*(Upsilon+1)+N_L+1;

% C1约束
A1 = [zeros(1,num_sensors*Upsilon+num_sensors*(Upsilon+1)+N_L), -1];
b1 = -0.9;

% C2约束, 同时在约束中把传感器的传输成功概率写进去
A2 = [];
P_sucess_sensor = zeros(num_sensors,Upsilon+1);
for i = 1:num_sensors
    A2 = [A2 eye(Upsilon)];
    P_sucess_sensor(i,:) = 1-(1-P_success(i)).^[0:1:Upsilon];
end
A2 = [A2 zeros(size(A2,1),n-size(A2,2))];
b2 = slot*ones(Upsilon,1);

% C3约束和C4约束
A3 = zeros(N_L,n); A3(:,end) = -ones(N_L,1);
A4 = zeros(N_L,n); A4(:,end) = ones(N_L,1);
b3 = (g_l-1)*ones(N_L,1);
b4 = -(g_l-1)*ones(N_L,1)+M;
ite = 1;
for j = 1:N_L
    for i = 1:g_l
        idx = L_mat(j,i);
        A3(j,num_sensors*Upsilon+(idx-1)*(Upsilon+1)+1:num_sensors*Upsilon+idx*(Upsilon+1)) = P_sucess_sensor(idx,:);
        A4(j,num_sensors*Upsilon+(idx-1)*(Upsilon+1)+1:num_sensors*Upsilon+idx*(Upsilon+1)) = -P_sucess_sensor(idx,:);
        A4(j,num_sensors*Upsilon+num_sensors*(Upsilon+1)+j) = M;
    end
end
% 总约束
A_total = [A1; A2; A3; A4];
b_total = [b1; b2; b3; b4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TII文章-面向可靠性的传感器调度
% 根据返修意见，增加强化学习、启发式算法、事件触发算法、凸集构建以及AoI方法的对比
% Peizhe Li
% 2025-07-01
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
lambda = 50;         % 热导率 W/(m·K)
delta_b = 0.05;      % 每个厚度单元的厚度 (m)
delta_l = 10;         % 每个长度单元长度 (m)
delta_t = 0.1;       % 时间步长 (s)
v = 0/60;               % 输送速度 (m/s)
% v=0
epsilon = 0.85;
sigma0 = 5.67e-8;
x_inf = 60 + 273.15;     % 恒定外界温度
T_total = 200/delta_t;   % 仿真总步数

% --------------------------
% 钢板空间离散结构
% --------------------------
nu = 2;                % 厚度方向分段数
tau_s = 5;            % 长度方向分段数
n = nu * tau_s;        % 状态变量维度

% --------------------------
% 传感器数量设定并初始化运行时间参数
% --------------------------
num_sensors = 10;
time_sub_avg = zeros(4,10);
time_opt_avg = zeros(1,10);
max_ite = 40;
time_sub_sum = [];
time_opt_sum = [];

x0 = 1180*ones(n,1);        % 初始温度为1180摄氏度
W = 0.1 * eye(n);       % 过程噪声协方差
V = 0.1 * eye(num_sensors);  % 观测噪声协方差

% --------------------------
% 优化问题参数
% --------------------------
zeta = 15;
Upsilon = 5;

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
SNR0 = 2000;     % SNR门限
% p_i = 10;     % 发射功率
% 理论概率
x = (N0 * SNR0 ) ./ (transmit_power * gains);
P_success = exp(-x); % 每个传感器的传输成功率

% --------------------------
% 设备温度模型参数
% 传感器可靠性参数，看xiaofan yu开源的代码确认一下
% --------------------------
aT = 0.9;
bT = 0.15;
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
% surface_choice = randi([1, 2], num_sensors, 1);  % 1 表示上表面，2 表示下表面
surface_choice = ones(num_sensors, 1);  % 1 表示上表面，2 表示下表面
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

[isObservable, ObsvMatrix, rankObsv, nStates] = checkObservability(F, C);

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

X_hist = []; gamma_hist = [];

%% 对subproblem数量进行循环
gamma1 = [];T_init_rec=[];
T_Simulation = 200;
for tq = 1:Upsilon:T_Simulation
    num_subproblem = 9;
    
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
    
    N_mod = mod(N_M,num_subproblem);
    L_mat = [L_mat;L_mat(1:num_subproblem-N_mod,:)];
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
    
    X_hist=[]; gamma_hist=[];
    
    % 开最大数量的并行
    parfor step = 1:num_subproblem
        Local_data = L_sub(start_point_L(step,1):end_point_L(step,1),:);    
        Local_A = A_sub(start_point_A(step,1):end_point_A(step,1),:);
        Local_b = b_sub(start_point_A(step,1):end_point_A(step,1),:);
        [X_value, f_value, time] = slove_X(Local_data, p, zeta, T_init, aT, bT, cT, Tm, num_sensors, Upsilon, Local_A, Local_b);
        X{step}= X_value; Obj_value{step} = f_value; 
    end
    
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

%% 执行AoI最优算法
AoI = zeros(num_sensors,1);
% AoI = [1:5 1:3 1:2];
% 生成一个调度序列
Success_sensor1 = rand(num_sensors, T_Simulation); 
for i = 1:num_sensors
    Success_sensor(i, Success_sensor1(i,:) > P_success(i)) = 0;
    Success_sensor(i, Success_sensor1(i,:) <= P_success(i)) = 1;
end

gamma_greedy_rec = [];
for tgreddy = 1:T_Simulation
    gamma_greedy = binvar(num_sensors,1);
    
    Objective = ones(1,num_sensors)*(diag(AoI)+eye(num_sensors))*(1-gamma_greedy)/num_sensors+zeta*p*sum(gamma_greedy)*delta_t;
    
    A1 = ones(1,num_sensors); b1 = slot;
    Constraints = A1*gamma_greedy <= b1;  
    
    options = sdpsettings('solver', 'gurobi', 'verbose', 1);
    result = optimize(Constraints, Objective, options);
    
    gamma_greedy_rec = [gamma_greedy_rec, value(gamma_greedy)];

    AoI = (AoI+1).*(1-and(gamma_greedy_rec(:,end),Success_sensor(:,tq)));
end

%% 执行能观性凸集对应算法
% 先有卡尔曼协方差的更新，然后用每个时刻的P的迭代去做优化
P = eye(size(F,1));
gamma_convex = zeros(num_sensors,T_Simulation);
% 在每一个时刻开始去做贪心
for tconvex = 1:T_Simulation
    % 协方差的一步预测
    P = F*P*F'+W;
    % 存储选择的传感器
    C_choice = [];
    % 初始化可用的传感器集合
    Sensor_set = 1:num_sensors;
    % 对每一个传感器都去做对应的调度
    for b = 1:slot
        % 初始化每个传感器选择后对应的P
        P_sensor = 10^3*ones(1,num_sensors);
        for i = Sensor_set
            C_choice_pre = [C_choice;C(i,:)];
            K = P*C_choice_pre'*inv(C_choice_pre*P*C_choice_pre'+eye(size(C_choice_pre,1)));
            P_sensor(1,i) = trace((eye(size(K,1))-K*C_choice_pre)*P);
        end
        [~,ind_min] = find(P_sensor==min(P_sensor));
        ind_min = ind_min(randi(numel(ind_min)));
        for i = 1:length(Sensor_set)
            if Sensor_set(i) == ind_min
                Sensor_set(i) = [];
                break;
            end
        end
        C_choice = [C_choice;C(ind_min,:)];
        gamma_convex(ind_min,tconvex) = 1;
    end
    K = P*C_choice'*inv(C_choice*P*C_choice'+eye(size(C_choice,1)));
    P = (eye(size(K,1))-K*C_choice)*P;
end

%% 执行事件触发算法
P0 = 25*eye(size(F,1));
x_hat_0 = zeros(n,1);
T_Simulation = 20000;
[X_true, X_hat, P4_trace, gamma_event] = simulate_event_triggered_multisensor_kf( ...
    F, C, W, V, x0, P0, x_hat_0, T_Simulation, 3, slot, P_success, 5);

[X_true, X_hat, P4_trace_1, gamma_event_1] = simulate_event_triggered_multisensor_kf( ...
    F, C, W, V, x0, P0, x_hat_0, T_Simulation, 10 ...
    , slot, P_success, 5);


%% 执行贪心算法
P0 = 25*eye(size(F,1));
x_hat_0 = zeros(n,1);
[X_true, X_hat, P5_trace, gamma_greedy] = simulate_greedy_slot_kf( ...
    F, C, W, V, x0, P0, x_hat_0, T_Simulation, 3, P_success, 45);
% plot(P5_trace)

T_Simulation = 200;
%% 强化学习算法
% -----------------------------
% Hyperparameters
% -----------------------------
episodes       = 4000;     % can increase if needed
horizon        = 200;
gamma          = 0.99;
alpha_lr       = 0.15;     % Q-learning step size
eps_start      = 1.0;
eps_end        = 0.1;
eps_decaySteps = 120000;   % in environment steps
AoI_cap        = 10;       % IMPORTANT: cap AoI for finite state space
print_every    = 200;

% -----------------------------
% Environment parameters
% -----------------------------
nSensors  = 10;
nSubsys   = 5;

% -----------------------------
% Build action set: all combinations of 5 sensors out of 10
% Each action is a 1x10 binary vector with exactly 5 ones
% -----------------------------
comb1 = nchoosek(1:nSensors, 2);     % (252 x 5)
comb2 = nchoosek(1:nSensors, 1);
nActions = size(comb1, 1) + size(comb2, 1);           % 252
A = zeros(nActions, nSensors);
for k = 1:size(comb1, 1)
    A(k, comb1(k,:)) = 1;
end
for k = size(comb1, 1)+1:nActions
    A(k, comb2(k-size(comb1, 1),:)) = 1;
end

fprintf('Action space: choose 3 of 10 -> nActions = %d\n', nActions);

% -----------------------------
% Q-table (sparse): containers.Map
% key: "a_b_c_d_e"  (capped AoI state)
% val: 1 x nActions Q-values
% -----------------------------
Q = containers.Map('KeyType','char','ValueType','any');

episode_return = zeros(episodes,1);
total_steps = 0;

% -----------------------------
% Training loop
% -----------------------------
for ep = 1:episodes
    % reset
    t = 0;
    aoi = zeros(1, nSubsys); % initial AoI = 0
    G = 0;
    done = false;

    while ~done
        eps = eps_schedule(total_steps, eps_start, eps_end, eps_decaySteps);

        % state key (capped)
        aoi_capd = min(max(round(aoi), 0), AoI_cap);
        sKey = stateKey(aoi_capd);

        % ensure Q row exists
        Qs = getQrow(Q, sKey, nActions);

        % epsilon-greedy over 252 actions
        if rand < eps
            aIdx = randi(nActions);
        else
            maxv = max(Qs);
            idxs = find(Qs == maxv);
            aIdx = idxs(randi(numel(idxs))); % random tie-break
        end

        % ---------- environment step ----------
        % activate 5 sensors according to action A(aIdx,:)
        active = A(aIdx,:); % 1x10 {0,1}

        % sample successes for active sensors
        success_mask = zeros(1, nSensors);
        actIdx = find(active == 1);
        u = rand(1, numel(actIdx));
        ok = (u < P_success(actIdx));
        ok = (u < 1);
        success_mask(actIdx) = ok;

        % update AoI per subsystem (2 sensors each)
        new_aoi = aoi;
        for i = 1:nSubsys
            s0 = 2*i - 1;
            s1 = 2*i;
            if success_mask(s0)==1 || success_mask(s1)==1
                new_aoi(i) = 0;
            else
                new_aoi(i) = new_aoi(i) + 1;
            end
        end
        aoi = new_aoi;

        % reward
        r = 1000 - sum( 1.01 .^ aoi ) - sum(active)*300;
 
        % time / done
        t = t + 1;
        done = (t >= horizon);

        % next state key
        aoi_next_capd = min(max(round(aoi), 0), AoI_cap);
        sNextKey = stateKey(aoi_next_capd);
        QsNext = getQrow(Q, sNextKey, nActions);

        % Q-learning update
        if done
            target = r;
        else
            target = r + gamma * max(QsNext);
        end
        Qs(aIdx) = Qs(aIdx) + alpha_lr * (target - Qs(aIdx));
        Q(sKey) = Qs; % write back

        G = G + r;
        total_steps = total_steps + 1;
    end

    episode_return(ep) = G;

    if mod(ep, print_every)==0
        avgR = mean(episode_return(ep-print_every+1:ep));
        fprintf('Episode %5d | avg return(last %d) = %.3f | eps = %.3f\n', ...
            ep, print_every, avgR, eps);
    end
end

% -----------------------------
% Evaluate greedy policy once (rollout) and export matrices
% -----------------------------
aoi = zeros(1,nSubsys);
aoi_traj = zeros(horizon, nSubsys);
action_traj = zeros(horizon, nSensors);
success_traj = zeros(horizon, nSensors);
reward_traj = zeros(horizon, 1);
action_index_traj = zeros(horizon, 1);

for k = 1:horizon
    aoi_capd = min(max(round(aoi), 0), AoI_cap);
    sKey = stateKey(aoi_capd);
    Qs = getQrow(Q, sKey, nActions);

    % greedy action (tie-break random)
    maxv = max(Qs);
    idxs = find(Qs == maxv);
    aIdx = idxs(randi(numel(idxs)));

    active = A(aIdx,:);
    actIdx = find(active == 1);

    success_mask = zeros(1, nSensors);
    u = rand(1, numel(actIdx));
    ok = (u < P_success(actIdx));
    success_mask(actIdx) = ok;

    new_aoi = aoi;
    for i = 1:nSubsys
        s0 = 2*i - 1; s1 = 2*i;
        if success_mask(s0)==1 || success_mask(s1)==1
            new_aoi(i) = 0;
        else
            new_aoi(i) = new_aoi(i) + 1;
        end
    end
    aoi = new_aoi;

    r = 1000 -sum( 1.01 .^ aoi ) - sum(active)*300;

    % record
    aoi_traj(k,:) = aoi;
    action_traj(k,:) = active;
    success_traj(k,:) = success_mask;
    reward_traj(k) = r;
    action_index_traj(k) = aIdx;
end

gamma_rl = action_traj';

%% 执行可靠性相关的算法
% 对subproblem数量进行循环
gamma8 = [];T_init_rec=[];
T_Simulation = 200;
for tq = 1:Upsilon:T_Simulation
    num_subproblem = 9;
    
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
    
    N_mod = mod(N_M,num_subproblem);
    L_mat = [L_mat;L_mat(1:num_subproblem-N_mod,:)];
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
    
    X_hist=[]; gamma_hist=[];
    
    % 开最大数量的并行
    parfor step = 1:num_subproblem
        Local_data = L_sub(start_point_L(step,1):end_point_L(step,1),:);    
        Local_A = A_sub(start_point_A(step,1):end_point_A(step,1),:);
        Local_b = b_sub(start_point_A(step,1):end_point_A(step,1),:);
        [X_value, f_value, time] = slove_X_Reliability(Local_data, p, zeta, T_init, aT, bT, cT, Tm, num_sensors, Upsilon, Local_A, Local_b);
        X{step}= X_value; Obj_value{step} = f_value; 
    end
    
    % 输出优化结果
    % 首先把 Obj_value 转换为数值数组
    Obj_numeric = cellfun(@(c) c, Obj_value);
    
    % 找到最小值的位置
    [~, idx_min] = min(Obj_numeric);
    
    % 输出对应的 X
    X_min = X{idx_min};    
    gamma8 = [gamma8 reshape(X{idx_min}(1:num_sensors*Upsilon,end), Upsilon, num_sensors)'];
    
    % 更新对应的初始温度
    gamma_now = reshape(X{idx_min}(1:num_sensors*Upsilon,end), Upsilon, num_sensors)';
    for i = 1:num_sensors
        T_init(i,1) = compute_temperature_after_upsilon(gamma_now(i,:), T_init(i,1), p, aT, bT, cT, Tm);
    end
    T_init_rec=[T_init_rec T_init];
end

%% 验证卡尔曼滤波本身是否收敛
P3 = eye(size(F,1));
P3_trace = zeros(1,5000);
for i = 1:7000
    P3 = F*P3*F'+W;
    K3 = P3*C'/(C*P3*C'+0.1*eye(size(C,1)));
    P3 = (eye(size(K3,1))-K3*C)*P3;
    P3_trace(1,i) = trace(P3);
end
figure(7)
% plot(P3_trace)
%% 画出相应的对比图
% 首先应该对比估计性能
% gamma1 = gamma;
T_Simulation = 200;
p = 0.05 + transmit_power;
power = p;
Gain_Simulation = 20000/T_Simulation;
gamma1 = kron(ones(1,Gain_Simulation),gamma1);
gamma2 = kron(ones(1,Gain_Simulation),gamma_greedy_rec);
gamma3 = kron(ones(1,Gain_Simulation),gamma_convex);
gamma4 = kron(ones(1,Gain_Simulation),gamma_event(:,1:T_Simulation));
gamma4_1 = kron(ones(1,Gain_Simulation),gamma_event_1(:,1:T_Simulation));
gamma5 = kron(ones(1,Gain_Simulation),gamma_greedy(:,1:T_Simulation));
gamma6 = kron(ones(1,Gain_Simulation),gamma_rl(:,1:T_Simulation));
gamma8 = kron(ones(1,Gain_Simulation),gamma8);
% clear gamma gamma_greedy_rec gamma_convex

P1 = 25*eye(size(F,1)); P1_trace = zeros(1,T_Simulation*Gain_Simulation);
P2 = 25*eye(size(F,1)); P2_trace = zeros(1,T_Simulation*Gain_Simulation);
P3 = 25*eye(size(F,1)); P3_trace = zeros(1,T_Simulation*Gain_Simulation);
P4 = 25*eye(size(F,1)); P4_trace = zeros(1,T_Simulation*Gain_Simulation);
P4_1 = 25*eye(size(F,1)); P4_trace_1 = zeros(1,T_Simulation*Gain_Simulation);
P5 = 25*eye(size(F,1)); P5_trace = zeros(1,T_Simulation*Gain_Simulation);
P6 = 25*eye(size(F,1)); P6_trace = zeros(1,T_Simulation*Gain_Simulation);
P8 = 25*eye(size(F,1)); P8_trace = zeros(1,T_Simulation*Gain_Simulation);

Tem1 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem1(:,1) = Tm*ones(num_sensors,1); %Tm*ones(num_sensors,1)
Tem2 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem2(:,1) = Tm*ones(num_sensors,1);
Tem3 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem3(:,1) = Tm*ones(num_sensors,1);
Tem4 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem4(:,1) = Tm*ones(num_sensors,1);
Tem4_1 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem4_1(:,1) = Tm*ones(num_sensors,1);
Tem5 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem5(:,1) = Tm*ones(num_sensors,1);
Tem6 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem6(:,1) = Tm*ones(num_sensors,1);
Tem8 = zeros(num_sensors,T_Simulation*Gain_Simulation); Tem8(:,1) = Tm*ones(num_sensors,1);

for i = 1:T_Simulation*Gain_Simulation
    % 把对应的观测矩阵拿出来
    idx = find(gamma1(:,i) ~= 0);
    C1 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C1 = [C1;C(idx(k),:)];
        end
    end

    idx = find(gamma2(:,i) ~= 0);
    C2 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C2 = [C2;C(idx(k),:)];
        end
    end

    idx = find(gamma3(:,i) ~= 0);
    C3 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C3 = [C3;C(idx(k),:)];
        end
    end

    idx = find(gamma4(:,i) ~= 0);
    C4 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C4 = [C4;C(idx(k),:)];
        end
    end

    idx = find(gamma4_1(:,i) ~= 0);
    C4_1 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C4_1 = [C4_1;C(idx(k),:)];
        end
    end

    idx = find(gamma5(:,i) ~= 0);
    C5 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C5 = [C5;C(idx(k),:)];
        end
    end

    idx = find(gamma6(:,i) ~= 0);
    C6 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C6 = [C6;C(idx(k),:)];
        end
    end

    idx = find(gamma8(:,i) ~= 0);
    C8 = [];
    transmit_index = rand(length(idx),1);
    for k = 1:length(idx)
        if transmit_index(k)<P_success(idx(k))
            C8 = [C8;C(idx(k),:)];
        end
    end

    % 协方差的一步预测
    P1 = F*P1*F'+W;
    P2 = F*P2*F'+W;
    P3 = F*P3*F'+W;
    P4 = F*P4*F'+W;
    P4_1 = F*P4_1*F'+W;
    P5 = F*P5*F'+W;
    P6 = F*P6*F'+W;
    P8 = F*P8*F'+W;

    
    % 协方差的更新
    if length(C1) > 0
        K1 = P1*C1'/(C1*P1*C1'+0.1*eye(size(C1,1)));
        P1 = (eye(size(K1,1))-K1*C1)*P1;
    end
    if length(C2) > 0
        K2 = P2*C2'/(C2*P2*C2'+0.1*eye(size(C2,1)));
        P2 = (eye(size(K2,1))-K2*C2)*P2;
    end
    if length(C3) > 0
        K3 = P3*C3'/(C3*P3*C3'+0.1*eye(size(C3,1)));
        P3 = (eye(size(K3,1))-K3*C3)*P3;
    end 
    if length(C4) > 0
        K4 = P4*C4'/(C4*P4*C4'+0.1*eye(size(C4,1)));
        P4 = (eye(size(K4,1))-K4*C4)*P4;
    end 
    if length(C4_1) > 0
        K4_1 = P4_1*C4_1'/(C4_1*P4_1*C4_1'+0.1*eye(size(C4_1,1)));
        P4_1 = (eye(size(K4_1,1))-K4_1*C4_1)*P4_1;
    end 
    if length(C5) > 0
        K5 = P5*C5'/(C5*P5*C5'+0.1*eye(size(C5,1)));
        P5 = (eye(size(K5,1))-K5*C5)*P5;
    end 
    if length(C6) > 0
        K6 = P6*C6'/(C6*P6*C6'+0.1*eye(size(C6,1)));
        P6 = (eye(size(K6,1))-K6*C6)*P6;
    end  
    if length(C8) > 0
        K8 = P8*C8'/(C8*P8*C8'+0.1*eye(size(C8,1)));
        P8 = (eye(size(K8,1))-K8*C8)*P8;
    end  

    % 协方差的迹存下来
    P1_trace(1,i) = trace(P1);
    P2_trace(1,i) = trace(P2);
    P3_trace(1,i) = trace(P3);
    P4_trace(1,i) = trace(P4);
    P4_trace_1(1,i) = trace(P4_1);
    P5_trace(1,i) = trace(P5);
    P6_trace(1,i) = trace(P6);
    P8_trace(1,i) = trace(P8);

    % 每个传感器温度的更新
    for j = 1:num_sensors
        Tem1(j,i+1) = compute_temperature_after_upsilon(gamma1(j,i), Tem1(j,i), power, aT, bT, cT, Tm);
        Tem2(j,i+1) = compute_temperature_after_upsilon(gamma2(j,i), Tem2(j,i), power, aT, bT, cT, Tm);
        Tem3(j,i+1) = compute_temperature_after_upsilon(gamma3(j,i), Tem3(j,i), power, aT, bT, cT, Tm);
        Tem4(j,i+1) = compute_temperature_after_upsilon(gamma4(j,i), Tem4(j,i), power, aT, bT, cT, Tm);
        Tem4_1(j,i+1) = compute_temperature_after_upsilon(gamma4_1(j,i), Tem4_1(j,i), power, aT, bT, cT, Tm);
        Tem5(j,i+1) = compute_temperature_after_upsilon(gamma5(j,i), Tem5(j,i), power, aT, bT, cT, Tm);
        Tem6(j,i+1) = compute_temperature_after_upsilon(gamma6(j,i), Tem6(j,i), power, aT, bT, cT, Tm);
        Tem8(j,i+1) = compute_temperature_after_upsilon(gamma8(j,i), Tem8(j,i), power, aT, bT, cT, Tm);
    end
end

%% 直接用平均温度去算可靠性增益
Tem1_avg = mean(Tem1,2); Tem2_avg = mean(Tem2,2); Tem3_avg = mean(Tem3,2); Tem4_avg_1 = mean(Tem4_1,2);
Tem4_avg = mean(Tem4,2); Tem5_avg = mean(Tem5,2); Tem6_avg = mean(Tem6,2); Tem8_avg = mean(Tem8,2); 
sensor_params = repmat(struct(), 1, num_sensors);
rng(130);
for i = 1:num_sensors
    sensor_params(i).delta = 2;
    sensor_params(i).kappa = 1.5e1 + 1*randn();
    sensor_params(i).J = 1.5e6 + 6*10^5*randn();
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
end
for i = 1:num_sensors
    p = sensor_params(i);
    T_avg(i) = compute_avg_temperature(p.aT, p.bT, p.cT, p.Tm, p.pi_t);
    % MTTF0(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, T_avg(i))/24/30/12;
    MTTF1(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem1_avg(i))/24/30/12;
    MTTF2(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem2_avg(i))/24/30/12;
    MTTF3(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem3_avg(i))/24/30/12;
    MTTF4(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem4_avg(i))/24/30/12;
    MTTF4_1(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem4_avg_1(i))/24/30/12;
    MTTF5(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem5_avg(i))/24/30/12;
    MTTF6(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem6_avg(i))/24/30/12;
    MTTF8(i) = estimate_mttf_from_Tavg(p.delta, p.kappa, p.J, p.Jc, p.Ea, p.B0, Tem8_avg(i))/24/30/12;
end

% 直接利用L_mat去计算系统的MTTF
[N_L,g_l] = size(L_mat);
MTTF_event = 10^2000*ones(8,N_L);
for i = 1:N_L 
    for j = 1:g_l
        MTTF_event(1,i) = min([MTTF_event(1,i),MTTF1(L_mat(i,j))]);
        MTTF_event(2,i) = min([MTTF_event(2,i),MTTF2(L_mat(i,j))]);
        MTTF_event(3,i) = min([MTTF_event(3,i),MTTF3(L_mat(i,j))]);
        MTTF_event(4,i) = min([MTTF_event(4,i),MTTF4(L_mat(i,j))]);
        MTTF_event(5,i) = min([MTTF_event(5,i),MTTF4_1(L_mat(i,j))]);
        MTTF_event(6,i) = min([MTTF_event(6,i),MTTF5(L_mat(i,j))]);
        MTTF_event(7,i) = min([MTTF_event(7,i),MTTF6(L_mat(i,j))]);
        MTTF_event(8,i) = min([MTTF_event(8,i),MTTF8(L_mat(i,j))]);
    end
end
MTTF_all1 = max(MTTF_event(1,:));
MTTF_all2 = max(MTTF_event(2,:));
MTTF_all3 = max(MTTF_event(3,:));
MTTF_all4 = max(MTTF_event(4,:));
MTTF_all4_1 = max(MTTF_event(5,:));
MTTF_all5 = max(MTTF_event(6,:));
MTTF_all6 = max(MTTF_event(7,:));
MTTF_all8 = max(MTTF_event(8,:));

% 还需要对比一下整体的能耗
Energy1 = sum(sum(gamma1))*transmit_power*delta_t;
Energy2 = sum(sum(gamma2))*transmit_power*delta_t;
Energy3 = sum(sum(gamma3))*transmit_power*delta_t;
Energy4 = sum(sum(gamma4))*transmit_power*delta_t;
Energy4_1 = sum(sum(gamma4_1))*transmit_power*delta_t;
Energy5 = sum(sum(gamma5))*transmit_power*delta_t;
Energy6 = sum(sum(gamma6))*transmit_power*delta_t;
Energy8 = sum(sum(gamma8))*transmit_power*delta_t;

% save('Comparison_Response_v5')
%% plot result
% 先对比估计效果
figure(2)
% gca.YScale = 'log';
hold on;

% plot([1:6000]*delta_t, P1_trace(1:6000), '-.', 'Color', [0.301, 0.745, 0.933], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'RoSS');
% 
% plot([1:6000]*delta_t, P2_trace(1:6000), '-.', 'Color', [0.466, 0.674, 0.188], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'Min-AoI');
% 
% plot([1:6000]*delta_t, P3_trace(1:6000), '-.', 'Color', [0.635, 0.078, 0.184], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'CSC');
% 
% plot([1:6000]*delta_t, P4_trace(1:6000), '-.', 'Color', [0.090, 0.235, 0.435], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'ET,$\sigma=3$');
% 
% plot([1:6000]*delta_t, P4_trace_1(1:6000), '-.', 'Color', [0.494, 0.698, 0.827], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'ET,$\sigma=10$');
% 
% plot([1:6000]*delta_t, P5_trace(1:6000), '-.', 'Color', [0.200, 0.600, 0.600], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'GT');
% 
% plot([1:6000]*delta_t, P6_trace(1:6000), '-.', 'Color', [0.850, 0.325, 0.098], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'RL');
% 
% plot([1:6000]*delta_t, P8_trace(1:6000), '-.', 'Color', [0.850, 0.325, 0.098], ...
%     'LineWidth', 1.8, 'MarkerSize', 8, 'DisplayName', 'MSR');

plot([1:6000]*delta_t, P1_trace(1:6000), '-', ...
    'Color', [0.000, 0.447, 0.741], ...   % 深蓝（主方法）
    'LineWidth', 2.0, 'DisplayName', 'RoSS');

plot([1:6000]*delta_t, P2_trace(1:6000), '--', ...
    'Color', [0.466, 0.674, 0.188], ...   % 绿色
    'LineWidth', 1.8, 'DisplayName', 'Min-AoI');

plot([1:6000]*delta_t, P3_trace(1:6000), ':', ...
    'Color', [0.635, 0.078, 0.184], ...   % 深红
    'LineWidth', 2.0, 'DisplayName', 'CSC');

plot([1:6000]*delta_t, P4_trace(1:6000), '-.', ...
    'Color', [0.301, 0.745, 0.933], ...   % 浅蓝（ET σ=1.5）
    'LineWidth', 1.8, 'DisplayName', 'ET,$\sigma=3$');

plot([1:6000]*delta_t, P4_trace_1(1:6000), '--', ...
    'Color', [0.090, 0.235, 0.435], ...   % 深蓝（ET σ=3）
    'LineWidth', 1.8, 'DisplayName', 'ET,$\sigma=10$');

plot([1:6000]*delta_t, P5_trace(1:6000), '-', ...
    'Color', [0.200, 0.600, 0.600], ...   % 青绿（GT）
    'LineWidth', 2.0, 'DisplayName', 'GT');

plot([1:6000]*delta_t, P6_trace(1:6000), '-', ...
    'Color', [0.850, 0.325, 0.098], ...   % 橙色（RL）
    'LineWidth', 2.0, 'DisplayName', 'RL');

plot([1:6000]*delta_t, P8_trace(1:6000), ':', ...
    'Color', [0.494, 0.184, 0.556], ...   % 紫色（MSR）
    'LineWidth', 2.0, 'DisplayName', 'MSR');

% set(gca, 'Interpreter','latex');
legend('Interpreter','latex');
set(gca, 'YScale', 'log');
% 坐标和标签
xlabel('Time (s)');
ylabel('tr(P) [K$^2$]','Interpreter','latex');
legend('Location', 'northwest');
set(gca,'FontName','Times New Roman','FontSize',16)
grid on;
box on;
% set(gca, 'FontSize', 15);

% ylim([0,20])

% 画图，整理成数组的形式，由于选择的weight不同，把对应的总cost的不要了
Bar_Data = [MTTF_all1/MTTF_all3;
            MTTF_all2/MTTF_all3;
            MTTF_all3/MTTF_all3;
            MTTF_all4/MTTF_all3;
            MTTF_all4_1/MTTF_all3;
            MTTF_all5/MTTF_all3;
            MTTF_all6/MTTF_all3;
            MTTF_all8/MTTF_all3;
            ];

Bar_Lab  = categorical({'RoSS','Min-AoI','CSC','ET,\sigma=3','ET,\sigma=10','GT','RL','MSR'});
Bar_Lab  = reordercats(Bar_Lab,{'RoSS','Min-AoI','CSC','ET,\sigma=3','ET,\sigma=10','GT','RL','MSR'});
figure(6)
set(gcf,'Position',[200 200 820 300]);  % [left bottom width height] 适当加宽
Cost_bar = bar(Bar_Lab,Bar_Data);
ylim([0.95, 1.06]);
Cost_bar(1).FaceColor = '#7fb8dd';
% Cost_bar(2).FaceColor = '#eba88b';
% Cost_bar(3).FaceColor = '#f5d78f';
% 显式替换 x 轴标签为 LaTeX 公式
% xticklabels({'RoSS','Min-AoI','CSC','ET,$\sigma=3$','ET,$\sigma=10$','GT','RL','MSR'})
xticklabels({
    'RoSS', ...
    'Min-AoI', ...
    'CSC', ...
    '$\begin{array}{c}\mathrm{ET}\\ \sigma=3\end{array}$', ...
    '$\begin{array}{c}\mathrm{ET}\\ \sigma=10\end{array}$', ...
    'GT', ...
    'RL', ...
    'MSR'
});
set(gca, 'TickLabelInterpreter', 'latex')
% ylabel('Cost')
set(gca,'FontName','Times New Roman','FontSize',12)
ylabel('Relative MTTF','FontName','Times New Roman','FontSize',16)
xlabel('Scheduling method','Interpreter','latex','FontName','Times New Roman','FontSize',16)
% legend('Estimation cost','Communication cost', 'Overall cost')
% legend('FontName','Times New Roman','FontSize',10.5,'LineWidth',1);

% 画图，能耗的对比
Bar_Data = [Energy1;
            Energy2;
            Energy3;
            Energy4;
            Energy4_1;
            Energy5;
            Energy6;
            Energy8
            ];

Bar_Lab  = categorical({'RoSS','Min-AoI','CSC','ET,\sigma=3','ET,\sigma=10','GT','RL','MSR'});
Bar_Lab  = reordercats(Bar_Lab,{'RoSS','Min-AoI','CSC','ET,\sigma=3','ET,\sigma=10','GT','RL','MSR'});
figure(7)
set(gcf,'Position',[200 200 820 300]);  % [left bottom width height] 适当加宽
Cost_bar = bar(Bar_Lab,Bar_Data);
% ylim([0.9, 1.2]);
Cost_bar(1).FaceColor = '#f5d78f';
% Cost_bar(2).FaceColor = '#eba88b';
% Cost_bar(3).FaceColor = '#f5d78f';

% 显式替换 x 轴标签为 LaTeX 公式
% xticklabels({'RoSS','Min-AoI','CSC','ET,$\sigma=3$','ET,$\sigma=10$','GT','RL','MSR'})
xticklabels({
    'RoSS', ...
    'Min-AoI', ...
    'CSC', ...
    '$\begin{array}{c}\mathrm{ET}\\ \sigma=3\end{array}$', ...
    '$\begin{array}{c}\mathrm{ET}\\ \sigma=10\end{array}$', ...
    'GT', ...
    'RL', ...
    'MSR'
});
set(gca, 'TickLabelInterpreter', 'latex')
% ylabel('Cost')
set(gca,'FontName','Times New Roman','FontSize',12)
ylabel('Energy consumption (mJ)','FontName','Times New Roman','FontSize',16)
xlabel('Scheduling method','Interpreter','latex','FontName','Times New Roman','FontSize',16)

% 传感器调度的间隔，一个传感器的调度，在这三个策略下的情况
% 计算当前点的x坐标宽度
gamma_all = [gamma1(5,:);gamma2(5,:);gamma3(5,:);gamma4(5,:);gamma4_1(5,:);gamma5(5,:);gamma6(5,:);gamma8(5,:)];
xWidth =1;
Height = 0.2;
Time_Int = 100;
figure(9)
for i = 1:Time_Int
    for j = 1:8
    % 绘制一个与x坐标宽度相同的小方块
        if gamma_all(j,i)==1
            patch([i, i, i+xWidth, i+xWidth], [j+Height, j-Height, j-Height, j+Height],[0.4660 0.6740 0.1880], 'FaceAlpha',.7);
            hold on;
        end
    end    
end
hold off;
box on
% legend({'$\alpha_1^\ast(t)$', '$\alpha_2^\ast(t)$', '$\alpha_3^\ast(t)$'}, 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold', 'NumColumns', 3);
axis([0 Time_Int,0.5,8+0.5])
yticks([1 2 3 4 5 6 7 8])
yticklabels({'RoSS','Min-AoI','CSC', 'ET,\sigma=3','ET,\sigma=10','GT','RL','MSR'})
xticks([0 20 40 60 80 100 120 140])
xticklabels({'0','2','4','6','8','10','12','14'})
xlabel('Time (s)','Interpreter','latex','fontsize',16,'FontWeight','bold')
% ylabel('Scheduling method','Interpreter','latex','fontsize',16,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',16)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 贪心的传感器调度算法
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_true, X_hat, traceP, gamma] = simulate_greedy_slot_kf( ...
    F, C, W, V, x0, P0, x_hat_0, T_Simulation, Slot, P_success, rng_seed)
%SIMULATE_GREEDY_SLOT_KF
% Greedy sensor selection (slot-limited) for remote KF:
% select Slot sensors that maximize expected trace reduction of covariance.
%
% Model:
%   x_{t+1} = F x_t + w_t, w_t ~ N(0,W)
%   y_t     = C x_t + v_t, v_t ~ N(0,V)
%   C: m x n, each row = one sensor (scalar measurement)
%
% Inputs:
%   F, C, W, V, x0, P0, x_hat_0, T_Simulation, Slot, P_success, rng_seed(optional)
%
% Outputs:
%   X_true : n x T_Simulation
%   X_hat  : n x T_Simulation
%   traceP : 1 x T_Simulation
%   gamma  : m x T_Simulation (1 if scheduled/triggered to transmit, else 0)

    if nargin >= 11 && ~isempty(rng_seed)
        rng(rng_seed);
    end

    [m, n] = size(C);
    assert(all(size(F)==[n,n]), 'F must be n x n');
    assert(all(size(W)==[n,n]), 'W must be n x n');
    assert(all(size(V)==[m,m]), 'V must be m x m');
    assert(all(size(x0)==[n,1]), 'x0 must be n x 1');
    assert(all(size(x_hat_0)==[n,1]), 'x_hat_0 must be n x 1');
    assert(all(size(P0)==[n,n]), 'P0 must be n x n');
    assert(length(P_success)==m, 'P_success must be length m');
    assert(Slot >= 0 && Slot <= m, 'Slot must be in [0,m]');

    % outputs
    X_true = zeros(n, T_Simulation);
    X_hat  = zeros(n, T_Simulation);
    traceP = zeros(1, T_Simulation);
    gamma  = zeros(m, T_Simulation);

    % init
    x = x0;
    x_hat = x_hat_0;
    P = P0;

    % noise sampling helpers
    LW = chol((W+W')/2 + 1e-12*eye(n), 'lower');
    LV = chol((V+V')/2 + 1e-12*eye(m), 'lower');

    for t = 1:T_Simulation
        % ---- True system evolve & measurement ----
        if t > 1
            w = LW * randn(n,1);
            x = F*x + w;
        end
        v = LV * randn(m,1);
        y = C*x + v;

        % ---- Remote KF: time update ----
        if t > 1
            x_hat = F*x_hat;
            P = F*P*F' + W;     % P^- (prior)
        end

        % ---- Greedy selection under Slot using expected trace reduction ----
        selected = greedy_select_sensors(P, C, V, P_success, Slot);

        % record scheduling decisions
        gamma(selected, t) = 1;

        % ---- Transmission + measurement updates (only if success) ----
        for kk = 1:length(selected)
            i = selected(kk);

            % success?
            % if rand() <= P_success(i)
            if rand() <= 1
                Ci = C(i,:);      % 1 x n
                Ri = V(i,i);      % scalar

                % scalar KF update
                nu = y(i) - Ci*x_hat;       % innovation
                S  = Ci*P*Ci' + Ri;         % scalar
                K  = (P*Ci')/S;             % n x 1

                I = eye(n);
                x_hat = x_hat + K*nu;
                % Joseph form (stable)
                P = (I - K*Ci)*P*(I - K*Ci)' + K*Ri*K';
            end
        end

        % ---- log ----
        X_true(:,t) = x;
        X_hat(:,t)  = x_hat;
        traceP(t)   = trace(P);
    end
end


function selected = greedy_select_sensors(Pprior, C, V, P_success, Slot)
%GREEDY_SELECT_SENSORS
% Greedy selection maximizing expected trace reduction:
% each step picks sensor i maximizing p_i*(tr(P)-tr(P_i_plus)),
% where P_i_plus is covariance after a scalar measurement update using sensor i.
%
% We update a "virtual" covariance Pvirt after each greedy pick to capture
% diminishing returns.

    [m, n] = size(C); %#ok<ASGLU>
    if Slot == 0
        selected = [];
        return;
    end

    remaining = true(m,1);
    selected  = zeros(Slot,1);

    Pvirt = Pprior;  % virtual covariance used for evaluating marginal gain

    for s = 1:Slot
        best_i = -1;
        best_gain = -inf;

        trP = trace(Pvirt);

        for i = 1:m
            if ~remaining(i), continue; end

            Ci = C(i,:);
            Ri = V(i,i);

            % Compute posterior covariance after using sensor i once
            S  = Ci*Pvirt*Ci' + Ri;      % scalar
            if S <= 0
                continue;
            end
            K  = (Pvirt*Ci')/S;          % n x 1
            I  = eye(size(Pvirt,1));
            Pplus = (I - K*Ci)*Pvirt*(I - K*Ci)' + K*Ri*K';  % Joseph form
            
            P_success1 = P_success;
            gain = (trP - trace(Pplus));      % expected trace reduction

            if gain > best_gain
                best_gain = gain;
                best_i = i;                 
            end
            if gain == best_gain
                sel = [best_i i];
                best_gain = gain;
                % best_i = sel(randi(numel(sel))); 
                if rand(1)<0.5
                    best_i = best_i;
                else
                    best_i = i;
                end
            end
        end

        if best_i < 0
            % no valid sensor found
            selected = selected(1:s-1);
            return;
        end

        selected(s) = best_i;
        remaining(best_i) = false;

        % Update virtual covariance as if this measurement is received
        Ci = C(best_i,:);
        Ri = V(best_i,best_i);
        S  = Ci*Pvirt*Ci' + Ri;
        K  = (Pvirt*Ci')/S;
        I  = eye(size(Pvirt,1));
        Pvirt = (I - K*Ci)*Pvirt*(I - K*Ci)' + K*Ri*K';
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 事件触发的传感器调度算法
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_true, X_hat, traceP, gamma] = simulate_event_triggered_multisensor_kf( ...
    F, C, W, V, x0, P0, x_hat_0, T_Simulation, Threshold, Slot, P_success, rng_seed)
%SIMULATE_EVENT_TRIGGERED_MULTISENSOR_KF
% Event-triggered multi-sensor remote state estimation with slot constraint.
%
% Inputs:
%   F            : n x n state transition matrix
%   C            : m x n measurement matrix (each row = one sensor, scalar measurement)
%   W            : n x n process noise covariance
%   V            : m x m measurement noise covariance (recommend diagonal if sensors independent)
%   x0           : n x 1 true initial state
%   P0           : n x n initial estimation error covariance at remote estimator
%   x_hat_0      : n x 1 initial state estimate at remote estimator
%   T_Simulation : simulation length (positive integer)
%   Threshold    : event-trigger threshold (scalar)
%   Slot         : max number of sensors allowed to attempt transmission per time step
%   P_success    : m x 1 vector, success probability for each sensor upon trigger
%   rng_seed     : (optional) integer for reproducibility
%
% Outputs:
%   X_true  : n x T_Simulation, true states
%   X_hat   : n x T_Simulation, remote estimator states
%   traceP  : 1 x T_Simulation, trace of estimation error covariance
%   gamma   : m x T_Simulation, trigger indicator (1=trigger, 0=no trigger)

    if nargin >= 12 && ~isempty(rng_seed)
        rng(rng_seed);
    end

    % Dimensions
    [m, n] = size(C);
    assert(all(size(F)==[n,n]), 'F must be n x n');
    assert(all(size(W)==[n,n]), 'W must be n x n');
    assert(all(size(V)==[m,m]), 'V must be m x m');
    assert(all(size(x0)==[n,1]), 'x0 must be n x 1');
    assert(all(size(x_hat_0)==[n,1]), 'x_hat_0 must be n x 1');
    assert(all(size(P0)==[n,n]), 'P0 must be n x n');
    assert(length(P_success)==m, 'P_success must be length m');

    % Pre-allocate outputs
    X_true = zeros(n, T_Simulation);
    X_hat  = zeros(n, T_Simulation);
    traceP = zeros(1, T_Simulation);
    gamma  = zeros(m, T_Simulation);

    % Initial conditions
    x = x0;
    x_hat = x_hat_0;
    P = P0;

    % Each sensor stores the last successfully delivered measurement to remote side
    % If never successfully sent before, we initialize it with its measurement at t=1
    y_last_success = nan(m,1);

    % Helper: sample Gaussian noise with covariance via Cholesky
    LW = chol((W+W')/2 + 1e-12*eye(n), 'lower');
    LV = chol((V+V')/2 + 1e-12*eye(m), 'lower');

    for t = 1:T_Simulation
        % ---- True system evolution & measurement ----
        if t == 1
            % x already x0
        else
            w = LW * randn(n,1);
            x = F * x + w;
        end

        v = LV * randn(m,1);
        y = C * x + v;          % m x 1, scalar per sensor

        % Initialize last_success at first step if needed
        if t == 1
            y_last_success = y; % treat initial as "known baseline" for triggering
        end

        % ---- Remote estimator: time update (prediction) ----
        if t == 1
            % keep initial x_hat_0, P0 as "time 1 prior" (or you can predict once here if desired)
        else
            x_hat = F * x_hat;
            P = F * P * F' + W;
        end

        % ---- Event trigger decision (based on last successful measurement) ----
        % delta = abs(y - y_last_success);          % m x 1
        delta = abs(y - C*x_hat);          % m x 1
        triggered = (delta > Threshold);          % m x 1 logical
        gamma(:, t) = double(triggered);

        trig_idx = find(triggered);
        % Apply slot constraint: at most Slot sensors can attempt transmission
        if ~isempty(trig_idx)
            if Slot < length(trig_idx)
                [~, order] = sort(delta(trig_idx), 'descend');
                trig_idx = trig_idx(order(1:Slot));
            end
        end

        % ---- Transmission success and KF measurement updates ----
        % We do scalar sequential updates for each successfully received sensor.
        % Only when transmission succeeds do we update y_last_success(i).
        for k = 1:length(trig_idx)
            i = trig_idx(k);

            % Bernoulli success
            if rand() <= P_success(i)
                % Measurement update with sensor i (scalar)
                Ci = C(i, :);          % 1 x n
                Ri = V(i, i);          % scalar (assume scalar noise per sensor)

                % Innovation
                yi = y(i);
                nu = yi - Ci * x_hat;  % scalar

                % Innovation covariance
                S  = Ci * P * Ci' + Ri;  % scalar

                % Kalman gain
                K  = (P * Ci') / S;       % n x 1

                % State and covariance update (Joseph form for numerical stability)
                I = eye(n);
                x_hat = x_hat + K * nu;
                P = (I - K*Ci) * P * (I - K*Ci)' + K*Ri*K';

                % Update last successful measurement baseline
                y_last_success(i) = yi;
            end
        end

        % ---- Log ----
        X_true(:, t) = x;
        X_hat(:, t)  = x_hat;
        traceP(t)    = trace(P);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 所提出算法相关函数定义
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

%% -----------------------------
% Local helper functions
% -----------------------------
function eps = eps_schedule(total_steps, eps_start, eps_end, eps_decaySteps)
    if total_steps >= eps_decaySteps
        eps = eps_end;
    else
        eps = eps_start + (eps_end - eps_start) * (total_steps / eps_decaySteps);
    end
end

function key = stateKey(aoi_vec)
    % Convert 1x5 AoI vector to a unique key "a_b_c_d_e"
    key = sprintf('%d_%d_%d_%d_%d', aoi_vec(1), aoi_vec(2), aoi_vec(3), aoi_vec(4), aoi_vec(5));
end

function Qrow = getQrow(Qmap, key, nActions)
    % Get Q(s,:) for key; init if new
    if isKey(Qmap, key)
        Qrow = Qmap(key);
    else
        Qrow = zeros(1, nActions);
        Qmap(key) = Qrow;
    end
end

%% 可靠性相关的算法求解，直接优化最大的传感器温度
%% initial value
function [X_value, f_value, time] = slove_X_Reliability(L_mat , p, zeta, T_init, a_T, b_T, c_T, T_m, num_sensors, Upsilon, A_total, b_total)
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
term1 = 0;
for i = 1:N
    temp_sum = 0;
    for k = 0:Upsilon-1
        for j = 0:k-1
            temp_sum = temp_sum + (a_T^(k-1-j)) * (b_T*p(i) * X((i-1)*Upsilon+1+j) + c_T + T_m);
        end
    end
    term1 = temp_sum / Upsilon;
    % term1 = max(term1,temp_sum);
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

% term1 = sum_T / N_L;
term1 = max(T_bar);

term2 = sum(sum(kron(p',ones(1,Upsilon))* X(1:num_sensors*Upsilon,:)));  % energy cost
% term2 = 0;  % energy cost
Objective = (term1 + 100 * term2);

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
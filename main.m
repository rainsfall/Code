% 基于Copula函数的风光功率联合场景生成
% 关键词：Copula；场景生成；风光出力相关性
clear; clc; close all;
%% 风光场景生成
%Copula;
%% 索引定义
i = 2;
t = 96; %时间点为96个点
dt = 60;    %时间间隔为15min


% %% 目标函数变量定义
% %   R_lm = sdpvar(1,96);中长期市场收益,决定变量为mu_nt_lm，P_ts_WPH_lm
% P_ts_WPH_lm = intvar(1,96); %中长期市场联合体中标电量
% 
% %   R_da = sdpvar(1,96);日前市场收益,决定变量为mu_ntw_da，P_wts_WPH_da
% mu_ntw_da = sdpvar(1,96);   %日前市场出清电价，也是下层功率平衡约束的拉格朗日乘子
% P_wts_WPH_da = intvar(1,96);    %日前市场联合体中标电量
% 
% %   R_rt = sdpvar(1,96);实时市场收益,决定变量为mu_ntw_rt，P_wts_WPH_rt
% mu_ntw_rt = sdpvar(1,96);   %实时市场出清电价，也是下层功率平衡约束的拉格朗日乘子
% P_wts_WPH_rt = intvar(1,96);    %实时市场联合体中标电量
% 
% %   R_TGC = sdpvar(1,96);绿证市场收益，决定变量为lambda_pt_TGC,Q_pt_TGC
% lambda_pt_TGC = intvar(1,1);    %当下市场绿证单价
% Q_pt_TGC = intvar(1,1); %当下市场绿证出售数量


%% 约束条件部分
%% 约束部分变量定义
V_it = sdpvar(2,96);
Q_it_pump = sdpvar(1,96);
Q_it = sdpvar(2,96);
Q_it_ele = sdpvar(2,96);
Q_it_dis = sdpvar(2,96);
P_itH = sdpvar(2,96);
P_itPump = sdpvar(1,96);

Constraints = [];
%% 水量平衡方程
mu = 700; % 均值
sigma = 1; % 标准差
I_it = max(normrnd(mu, sigma, [1, 96]), 0);

for t = 1:95
    % 上级水库的水量平衡约束条件
    Constraints = [Constraints, V_it(1,t+1) - V_it(1,t) - (I_it(1,t) + Q_it_pump(t) - Q_it_ele(1,t) - Q_it_dis(1,t)) * dt *60  == 0];
    Constraints = [Constraints, V_it(2,t+1) - V_it(2,t) - (Q_it(1,t) - Q_it_pump(t) - Q_it_ele(2,t) - Q_it_dis(2,t)) * dt *60  == 0];
end

%% 发电流量约束
% 下泄流量构成约束
Constraints = [Constraints,Q_it_ele + Q_it_dis == Q_it];
% 水位库容关系线性化
% Z_1tl_min = 2530;   %龙羊峡水电站死水位
% Z_1tl_max = 2600;   %龙羊峡水电站正常蓄水位
% Z_2tl_min = 2440;   %拉西瓦水电站死水位
% Z_2tl_max = 2452;   %拉西瓦水电站正常蓄水位

intervals1 = [2530, 2540, 2550, 2560, 2570, 2580, 2590, 2600]; % 龙羊峡水位
ref_values1 = [52.85, 71.90, 93.81, 117.62, 145.24, 175.71, 210.24, 247.14]*1e8; % 龙羊峡库容
intervals2 = [2440, 2445, 2449, 2452]; % 拉西瓦水位
ref_values2 = [8.55, 9.1, 9.5, 10]*1e8; % 拉西瓦库容

Z_it_upper = sdpvar(2, 96);

% 龙羊峡水位库容关系线性化
V_itl1 = sdpvar(1, 96, length(intervals1));
phi_itl1 = binvar(1, 96, length(intervals1));

for t = 1:96
    Constraints = [Constraints, sum(phi_itl1(1, t, :)) == 1];
    Constraints = [Constraints, sum(V_itl1(1, t, :)) == V_it(1, t)];
    
    local_Z = 0;
    
    for l = 1:length(intervals1)-1
        Constraints = [Constraints, implies(ref_values1(l) <= V_itl1(1, t, l) <= ref_values1(l+1), phi_itl1(1, t, l))];
        Constraints = [Constraints, phi_itl1(1, t, l) * ref_values1(l) <= V_itl1(1, t, l) <= phi_itl1(1, t, l) * ref_values1(l+1)];
        
        local_Z = local_Z + (phi_itl1(1, t, l) * intervals1(l) + (intervals1(l+1) - intervals1(l)) / (ref_values1(l+1) - ref_values1(l)) * (V_itl1(1, t, l) - phi_itl1(1, t, l) * ref_values1(l)));
    end
    
    Constraints = [Constraints, Z_it_upper(1, t) == local_Z];
end

% 拉西瓦水位库容关系线性化
V_itl2 = sdpvar(1, 96, length(intervals2));
phi_itl2 = binvar(1, 96, length(intervals2));
for t = 1:96
    local_constraints = [];
    local_Z = 0; % 局部变量存储 Z 的计算结果
    for l = 1:length(intervals2)-1
        local_constraints = [local_constraints, implies(ref_values2(l) <= V_itl2(1, t, l) <= ref_values2(l+1),phi_itl2(1, t, l))];
        local_constraints = [local_constraints,phi_itl2(1, t, l)*ref_values2(l) <= V_itl2(1, t, l) <= phi_itl2(1, t, l)*ref_values2(l+1)]
        local_Z = local_Z + (phi_itl2(1, t, l) * intervals2(l) + (intervals2(l+1) - intervals2(l)) / (ref_values2(l+1) - ref_values2(l)) * (V_itl2(1, t, l) - phi_itl2(1, t, l) * ref_values2(l)));
    end
    local_constraints = [local_constraints, sum(phi_itl2(1, t, :)) == 1];
    local_constraints = [local_constraints, sum(V_itl2(1, t, :)) == V_it(2, t)];
    local_constraints = [local_constraints, Z_it_upper(2, t) == local_Z]; % 将局部变量的值赋给全局变量
Constraints = [Constraints, local_constraints];
end

%尾水位——泄流量关系线性化
intervals3 = [2451,2452,2453,2454,2455]; % 龙羊峡尾水位
ref_values3 = [147.6215585,393.4826246,636.522124,879.1768642,1120.264067]; % 龙羊峡泄流量
intervals4 = [2236, 2237, 2238, 2239, 2240, 2241]; % 拉西瓦尾水位
ref_values4 = [791.7494505, 1172.953846, 1559.934066, 2007.23956,2460.962637,2917.252747]; % 拉西瓦泄流量

Z_it_down = sdpvar(2, 96);

%龙羊峡尾水位——泄流量关系线性化
Q_itm3 = sdpvar(1, 96, length(intervals3));
phi_itm3 = binvar(1, 96, length(intervals3));
for t = 1:96
    local_constraints = [];
    local_Z = 0; % 局部变量存储 Z 的计算结果
    for m = 1:length(intervals3)-1
        local_constraints = [local_constraints, implies(ref_values3(m) <= Q_itm3(1, t, m) <= ref_values3(m+1),phi_itm3(1, t, m))];
        local_constraints = [local_constraints,phi_itm3(1, t, m)*ref_values3(m) <= Q_itm3(1, t, m) <= phi_itm3(1, t, m)*ref_values3(m+1)]
        local_Z = local_Z + (phi_itm3(1, t, m) * intervals3(m) + (intervals3(m+1) - intervals3(m)) / (ref_values3(m+1) - ref_values3(m)) * (Q_itm3(1, t, m) - phi_itm3(1, t, m) * ref_values3(m)));
    end
    local_constraints = [local_constraints, sum(phi_itm3(1, t, :)) == 1];
    local_constraints = [local_constraints, sum(Q_itm3(1, t, :)) == Q_it(1, t)];
    local_constraints = [local_constraints, Z_it_down(1, t) == local_Z]; % 将局部变量的值赋给全局变量
Constraints = [Constraints, local_constraints];
end
%拉西瓦尾水位——泄流量关系线性化
Q_itm4 = sdpvar(1, 96, length(intervals4));
phi_itm4 = binvar(1, 96, length(intervals4));
for t = 1:96
    local_constraints = [];
    local_Z = 0; % 局部变量存储 Z 的计算结果
    for m = 1:length(intervals4)-1
        local_constraints = [local_constraints, implies(ref_values4(m) <= Q_itm4(1, t, m) <= ref_values4(m+1),phi_itm4(1, t, m))];
        local_constraints = [local_constraints,phi_itm4(1, t, m)*ref_values4(m) <= Q_itm4(1, t, m) <= phi_itm4(1, t, m)*ref_values4(m+1)]
        local_Z = local_Z + (phi_itm4(1, t, m) * intervals4(m) + (intervals4(m+1) - intervals4(m)) / (ref_values4(m+1) - ref_values4(m)) * (Q_itm4(1, t, m) - phi_itm4(1, t, m) * ref_values4(m)));
    end
    local_constraints = [local_constraints, sum(phi_itm4(1, t, :)) == 1];
    local_constraints = [local_constraints, sum(Q_itm4(1, t, :)) == Q_it(2, t)];
    local_constraints = [local_constraints, Z_it_down(2, t) == local_Z]; % 将局部变量的值赋给全局变量
Constraints = [Constraints, local_constraints];
end

% 水电站功率约束
Constraints = [Constraints,...
    P_itH(1,:) == 4*9.81*0.5*(Z_it_upper(1,:)-Z_it_down(1,:)).*Q_it_ele(1,:),...
    P_itH(2,:) == 7*9.81*0.5*(Z_it_upper(2,:)-Z_it_down(2,:)).*Q_it_ele(2,:)];
% 泵站耗电量约束
Constraints = [Constraints,...
    P_itPump == 4*9.81*0.5*(Z_it_upper(1,:)-Z_it_upper(2,:)).*Q_it_pump];

%% 水电站运行边界条件
% 水位约束
Constraints = [Constraints, 2530 <= Z_it_upper(1,:) <= 2600, 2440<=Z_it_upper(2,:)<=2452];
Constraints = [Constraints, 2451<=Z_it_down(1,:)<=2455, 2236<=Z_it_down(2,:)<=2241];
% 下泄流量约束
Constraints = [Constraints, 0<= Q_it(1,:) <= 1200, 0<= Q_it_ele(1,:) <= 1192, 0<=Q_it_dis(1,:)<=1000 ];
Constraints = [Constraints, 0<= Q_it(2,:) <= 2660, 0<= Q_it_ele(2,:) <= 2660, 0<=Q_it_dis(2,:)<=2000 ];
Constraints = [Constraints,0<= Q_itm3 <= 1192,0<= Q_itm4 <= 2660];
% 机组出力约束
Constraints = [Constraints, 589.8*1e3<= P_itH(1,:)<= 1280*1e3,990*1e3<=P_itH(2,:)<=4200*1e3];
% 泵站出力约束
Constraints = [Constraints,0<=Q_it_pump<=172.5*4];
% 库容约束
Constraints = [Constraints, 52.85*1e8<=V_it(1,:)<=247.14*1e8, 8.55*1e8<=V_it(2,:)<=10*1e8];
Constraints = [Constraints, 0*1e8<=V_itl1<=247.14*1e8, 0<=V_itl2<=10*1e8];
% 初值约束
Constraints = [Constraints,V_it(1,1) == 60*1e8,V_it(2,1) == 9.1*1e8]; %库容初值
Constraints = [Constraints,Q_it(1,1) == 500,Q_it(2,1) == 1559.93]; %泄流量初值

% 绝对值问题
load = normrnd(4*1e6, 1e6, [1, 96]);
Z = sdpvar(1,96);
Constraints = [Constraints,Z >= P_itH(1,:)+P_itH(2,:)-load(1,:),...
    Z >= load(1,:)-P_itH(1,:)-P_itH(2,:)];


ops = sdpsettings('solver', 'Gurobi+', 'verbose', 2, 'debug', 1);
ops.gurobi.IntegralityFocus = 1;
ops.gurobi.gurobi.NonConvex = 2;
obj = sum(Z) ;
result = optimize(Constraints,obj,ops);
P_itH_value = value(P_itH);
Q_it_value = value(Q_it);
Q_it_dis_value = value(Q_it_dis);
Q_it_ele_value = value(Q_it_ele);
Q_it_pump_value = value(Q_it_pump);
V_it_value = value(V_it);
P_itPump_value = value(P_itPump);
Z_it_upper_value = value(Z_it_upper);
obj_value = value(obj);

%机组发电的数据需要进一步验证。因为一天之内上下水位变化量实在太小了，按理来说不应该这么小的，注意观察
if result.problem ~= 0
    [model,recoverymodel,diagnostic,internalmodel] = export(Constraints, obj,sdpsettings('solver','GUROBI+'));
    iis = gurobi_iis(model);
    gurobi_write(model, 'TestModel.lp');
    find(iis.Arows);
end
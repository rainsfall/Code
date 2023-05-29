% 基于Copula函数的风光功率联合场景生成
% 关键词：Copula；场景生成；风光出力相关性
%clear; clc; close all;
%% 风光场景生成
%Copula;
%% 索引定义
i = 2;
t = 96; %时间点为96个点
dt = 15;    %时间间隔为15min


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
mu = 25; % 均值
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

Z_it_upper = sdpvar(2,96);
Z_it_down = sdpvar(2,96);
% 上下水位关系约束
polyfit_relation
Constraints = [Constraints,...
    Z_it_upper(1,:) == f1(V_it(1,:)*1e-8),...
    Z_it_upper(2,:) == f2(V_it(2,:)*1e-8),...
    Z_it_down(1,:) == f3(Q_it(1,:)*1e-8),...
    Z_it_down(2,:) == f4(Q_it(2,:)*1e-8)];
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
Constraints = [Constraints, 0<= Q_it(1,:) <= 1500, 0<= Q_it_ele(1,:) <= 1192];
Constraints = [Constraints, 0<= Q_it(2,:) <= 3000, 0<= Q_it_ele(2,:) <= 2660];

% 机组出力约束
Constraints = [Constraints, 589.8<= P_itH(1,:)*1e-3<= 1280,990<=P_itH(2,:)*1e-3<=4200];
% 泵站出力约束
Constraints = [Constraints,0<=Q_it_pump<=172.5*4];
% 库容约束
Constraints = [Constraints, 0<=V_it(1,:)<=1, 0<=V_it(2,:)<=1];
% 初值约束
Constraints = [Constraints,V_it(1,1) == 0.2,V_it(2,1) == 0.6]; %库容初值
Constraints = [Constraints,Q_it(1,1) == 393,Q_it(2,1) == 1559.93]; %泄流量初值

ops = sdpsettings('solver', 'Gurobi+', 'verbose', 2, 'debug', 1);
ops.gurobi.IntegralityFocus = 1;
ops.gurobi.gurobi.NonConvex = 2;
obj = sum(P_itH(1,:)+P_itH(2,:) - P_itPump) ;
result = optimize(Constraints,-obj,ops);
P_itH_value = value(P_itH);
Q_it_value = value(Q_it);
Q_it_dis_value = value(Q_it_dis);
Q_it_ele_value = value(Q_it_ele);
Q_it_pump_value = value(Q_it_pump);
V_it_value = value(V_it);
P_itPump_value = value(P_itPump);
Z_it_upper_value = value(Z_it_upper);

%机组发电的数据需要进一步验证。因为一天之内上下水位变化量实在太小了，按理来说不应该这么小的，注意观察
if result.problem ~= 0
    [model,recoverymodel,diagnostic,internalmodel] = export(Constraints, obj,sdpsettings('solver','GUROBI+'));
    iis = gurobi_iis(model);
    gurobi_write(model, 'TestModel.lp');
    find(iis.Arows);
end
% 基于Copula函数的风光功率联合场景生成
% 关键词：Copula；场景生成；风光出力相关性
clear; clc; close all;
%% 导入数据与预处理
solardata =  xlsread('数据-光伏.xlsx');
winddata = xlsread('数据-风功率.xlsx');
winddata = winddata(2:end, :);
solardata = solardata(2:end, :);

scenarionum = 500;  % 初始场景数目，可修改
num_cluster = 5;     % 要削减到的场景数目，可修改
ntime = 24;  % 24小时

% X和Y分别存储风和光的24个时刻历史观测数据
X = []; Y = [];
for t = 1 : ntime
    X{t} = winddata(:, t);
    Y{t} = solardata(:, t);
end


%% Copula拟合
% Frank-Copula 函数可以同时考虑变量的非负与负相关的关系
% 故采用 Frank-Copula 函数分别对24个时刻进行拟合

for i = 1 : ntime
    U = ksdensity(X{i}, 'function', 'cdf'); % 核密度估计
    V = ksdensity(Y{i}, 'function', 'cdf');
    
    alpha = copulafit('frank', [U(:) V(:)]); % 拟合出的参数
    copulaparams.alpha = alpha;
    copulaparams.numParams = 1;
    copModels(i) = copulaparams;       
end
%% 绘制二元Frank-Copula的密度函数和分布函数图
[Udata, Vdata] = meshgrid(linspace(0,1,31));  % 为绘图需要，产生新的网格数据
Ccdf_Frank = copulacdf('Frank', [Udata(:), Vdata(:)], copModels(12).alpha);


%% 采样
Data = cell(1, ntime); 

for j = 1: ntime
    Data{j} = copularnd('Frank', copModels(i).alpha, scenarionum); % 从Frank Copula中采样
end

%　逆变换，转换为实际场景
w = cell(1, ntime);
for k = 1 : ntime
    [fwind, xwind] = ecdf(X{k});
    funwind = @(pdwind)spline(fwind, xwind, pdwind);  % 样条插值
    
    [fsolar, xsolar] = ecdf(Y{k});
    funsolar = @(pdsolar)spline(fsolar, xsolar, pdsolar);
    
    for i = 1 : 2
        for j = 1 : scenarionum
            if i == 1
                w{k}(j, i) = funwind(Data{k}(j, i));  % 逆变换
            else
                w{k}(j,i) = funsolar(Data{k}(j, i));
            end
        end
    end
end
%% 场景削减
wind = [];
solar = [];

for i = 1 : ntime
    wind = [wind, w{i}(:, 1)];
    solar = [solar, w{i}(:, 2)];
end

% K-means削减
opts = statset('Display', 'final');
[idx, C] = kmeans([wind, solar], num_cluster, 'Distance', 'cityblock', 'Replicates', 20, 'Options',opts);

C_wind = C(:, 1:24);     % 最终的风电场景
C_solar = C(:, 25:end);  % 最终的光伏场景
%% 概率计算
for i = 1 : num_cluster
    p(i) =  length(find(idx == i)) / length(idx);
end
disp(['各场景概率为: ', num2str(p)]);

P_wd2=zeros(5,24);
for j=1:5
    for i=1:24
         P_wd2(j,i)=P_wd2(j,i)+p(j)*C_wind(j,i);
    end
end
P_WD=zeros(1,24);
for i=1:24
    for j=1:5
         P_WD(1,i)=P_WD(1,i)+P_wd2(j,i);
    end
end

P_pv2=zeros(5,24);
for j=1:5
    for i=1:24
         P_pv2(j,i)=P_pv2(j,i)+p(j)*C_solar(j,i);
    end
end

P_PV=zeros(1,24);
for i=1:24
    for j=1:5
         P_PV(1,i)=P_PV(1,i)+P_pv2(j,i);
    end
end
% 原始场景集画图
figure(1);
[ss,gg]=meshgrid(1:scenarionum,1:24);
plot3(ss,gg,wind, 'linewidth', 1);
title(['考虑相关性生成的风电出力', num2str(scenarionum), '个场景'])
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(2);
plot3(ss,gg,solar,'linewidth', 1);
title(['考虑相关性生成的光伏出力', num2str(scenarionum), '个场景'])
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');
% 画图
figure(3);
[ss1,gg1]=meshgrid(1:num_cluster,1:24);
plot3(ss1,gg1,C_wind, 'linewidth', 1);
title(['削减后风电出力', num2str(num_cluster), '个场景'])
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(4);
plot3(ss1,gg1,C_solar, 'linewidth', 1);
title(['削减后光伏出力', num2str(num_cluster), '个场景'])
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');
figure(5); 
bar(p); 
xlabel('场景数'); 
ylabel('概率');
title('各场景下概率');
figure(6)
plot(P_WD,'LineWidth',1.5); 
xlabel('时间/t'); 
ylabel('功率/kW');
title('	风电不确定性出力');
grid on
figure(7)
plot(P_PV,'LineWidth',1.5); 
xlabel('时间/t'); 
ylabel('功率/kW');
title('	光电不确定性出力');
grid on
figure(8);  
surf(Udata, Vdata, reshape(Ccdf_Frank,size(Udata)));  % 绘制二元Frank-Copula分布函数图
xlabel('风电出力'); 
ylabel('光伏出力'); 
zlabel('C(u,v)');  
title('二元Frank-Copula分布函数图')

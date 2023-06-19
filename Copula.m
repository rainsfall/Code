%% 导入数据与预处理
solardata =  csvread('光伏原出力数据.csv');
winddata = csvread('风电原出力数据.csv');
solardata_season = zeros(90,96,4);
winddata_season = zeros(90,96,4);
for i = 1:4
    solardata_temp = zeros(8640,1);
    winddata_temp = zeros(8640,1);
    solardata_temp(:,:) = solardata((i-1)*90*96+1:i*90*96,1);
    solardata_season(:,:,i) = reshape(solardata_temp,[96,90,1]).';
    winddata_temp(:,:) = winddata((i-1)*90*96+1:i*90*96,1);
    winddata_season(:,:,i) = reshape(winddata_temp,[96,90,1]).';
end
clear winddata_temp solardata_temp;
solardata_spring = solardata_season(:,:,1);
winddata_spring = winddata_season(:,:,1);
solardata_summer = solardata_season(:,:,2);
winddata_summer = winddata_season(:,:,2);
solardata_autumn = solardata_season(:,:,3);
winddata_autumn = winddata_season(:,:,3);
solardata_winter = solardata_season(:,:,4);
winddata_winter = winddata_season(:,:,4);
clear winddata_season solardata_season winddata solardata

%% 风光场景生成，分为四季
[data,time]=meshgrid(1:5,1:96);

[solardata_spring,winddata_spring] = Copula_function(solardata_spring,winddata_spring);
figure(1);
plot3(data,time,winddata_spring, 'linewidth', 1);
title('春天风电出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(2);
plot3(data,time,solardata_spring, 'linewidth', 1);
title('春天光伏出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');

[solardata_summer,winddata_summer] = Copula_function(solardata_summer,winddata_summer);
figure(3);
plot3(data,time,winddata_summer, 'linewidth', 1);
title('夏天风电出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(4);
plot3(data,time,solardata_summer, 'linewidth', 1);
title('夏天光伏出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');

[solardata_autumn,winddata_autumn] = Copula_function(solardata_autumn,winddata_autumn);
figure(5);
plot3(data,time,winddata_autumn, 'linewidth', 1);
title('秋天风电出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(6);
plot3(data,time,solardata_autumn, 'linewidth', 1);
title('秋天光伏出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');


[solardata_winter,winddata_winter] = Copula_function(solardata_winter,winddata_winter);
figure(7);
plot3(data,time,winddata_winter, 'linewidth', 1);
title('冬天风电出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('风电出力值');
figure(8);
plot3(data,time,solardata_winter, 'linewidth', 1);
title('冬天光伏出力5个代表场景')
xlabel('场景'); ylabel('时刻');zlabel('光伏出力值');

solardata_min_spring = min(solardata_spring);
solardata_max_spring = max(solardata_spring);
winddata_min_spring = min(winddata_spring);
winddata_max_spring = max(winddata_spring); 

solardata_min_summer = min(solardata_summer);
solardata_max_summer = max(solardata_summer);
winddata_min_summer = min(winddata_summer);
winddata_max_summer = max(winddata_summer);

solardata_min_autumn = min(solardata_autumn);
solardata_max_autumn = max(solardata_autumn);
winddata_min_autumn = min(winddata_autumn);
winddata_max_autumn = max(winddata_autumn);

solardata_min_winter = min(solardata_winter);
solardata_max_winter = max(solardata_winter);
winddata_min_winter= min(winddata_winter);
winddata_max_winter = max(winddata_winter);
%% 函数定义
function [C_solar,C_wind] = Copula_function(solardata,winddata)
    %数据扩大化
    solardata = solardata.*10;
    winddata = winddata.*10; 
    
    scenarionum = 500;  % 初始场景数目，可修改
    num_cluster = 5;     % 要削减到的场景数目，可修改
    ntime = 96;  % 96个点
    
    % X和Y分别存储风和光的96个时刻历史观测数据
    X = []; Y = [];
    for t = 1 : ntime
        X{t} = winddata(:, t);
        Y{t} = solardata(:, t);
    end
    %% Copula拟合
    % Frank-Copula 函数可以同时考虑变量的非负与负相关的关系
    % 故采用 Frank-Copula 函数分别对96个时刻进行拟合
    
    for i = 1 : ntime
        U = ksdensity(X{i}, 'function', 'cdf'); % 核密度估计
        V = ksdensity(Y{i}, 'function', 'cdf');
        
        alpha = copulafit('Frank', [U(:) V(:)]); % 拟合出的参数
        copulaparams.alpha = alpha;
        copulaparams.numParams = 1;
        copModels(i) = copulaparams;       
    end
    %% 绘制二元Frank-Copula的密度函数和分布函数图
    [Udata, Vdata] = meshgrid(linspace(0,1,31));  % 为绘图需要，产生新的网格数据
    Ccdf_Frank = copulacdf('Frank', [Udata(:), Vdata(:)], copModels(48).alpha);
    
    
    %% 采样
    Data = cell(1, ntime); 
    
    for j = 1: ntime
        Data{j} = copularnd('Frank', copModels(j).alpha, scenarionum); % 从Frank Copula中采样
    end
    
    %　逆变换，转换为实际场景
    w = cell(1, ntime);
    for k = 1 : ntime
        [fwind, xwind] = ecdf(X{k});
        funwind = @(pdwind)spline(fwind, xwind, pdwind);  % 样条插值
        
        [fsolar, xsolar] = ecdf(Y{k});
        funsolar = @(pdsolar)spline(fsolar, xsolar, pdsolar);
        
        for j = 1 : scenarionum
            u = Data{k}(j, 1);
            v = Data{k}(j, 2);
            w{k}(j, 1) = max(funwind(u), 0);    % 确保功率值非负
            w{k}(j, 2) = max(funsolar(v), 0);
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
    
    C_wind = C(:, 1:96);     % 最终的风电场景
    C_solar = C(:, 97:end);  % 最终的光伏场景
    %% 概率计算
    for i = 1 : num_cluster
        p(i) =  length(find(idx == i)) / length(idx);
    end
    disp(['各场景概率为: ', num2str(p)]);
    
    P_wd2=zeros(5,96);
    for j=1:5
        for i=1:96
             P_wd2(j,i)=P_wd2(j,i)+p(j)*C_wind(j,i);
        end
    end
    P_WD=zeros(1,96);
    for i=1:96
        for j=1:5
             P_WD(1,i)=P_WD(1,i)+P_wd2(j,i);
        end
    end
    
    P_pv2=zeros(5,96);
    for j=1:5
        for i=1:96
             P_pv2(j,i)=P_pv2(j,i)+p(j)*C_solar(j,i);
        end
    end
    
    P_PV=zeros(1,96);
    for i=1:96
        for j=1:5
             P_PV(1,i)=P_PV(1,i)+P_pv2(j,i);
        end
    end
    clear i time;
end

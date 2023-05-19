dataquanzhou = readmatrix("2016-2018负荷天气data-quanzhou.csv");
dataquanzhou = dataquanzhou(:,2);
load_data = reshape(dataquanzhou,[1096,96]);
load_data = load_data(1:92,:);
num_scenarios = 5; % 更改为5
num_samples = 1000;

% 拉丁超立方采样生成采样点
X = lhsdesign(num_samples, size(load_data, 2), 'criterion', 'maximin');
X = repmat(min(load_data), num_samples, 1) + X.*repmat(max(load_data) - min(load_data), num_samples, 1);

% 运用kmeans进行聚类，得到聚类中心
[idx, C] = kmeans(X, num_scenarios);

% 生成场景数据
scenario_data = zeros(num_scenarios, size(load_data, 2));
for i = 1:num_scenarios
    scenario_data(i,:) = C(idx(i), :);
end

% 计算每个场景的概率
counts = histcounts(idx, num_scenarios);
probabilities = counts / num_samples;

% 可视化场景数据
figure(1);
[datas,times]=meshgrid(1:num_scenarios,1:96);
plot3(datas,times,scenario_data(:,:));
title(['负荷功率', num2str(num_scenarios), '个场景'])
xlabel('场景'); ylabel('时刻');zlabel('负荷功率');

% 可视化某一天的数据
figure(2);
hold on;
for i = 1:num_scenarios
    plot(scenario_data(i,1:96), 'LineWidth', probabilities(i)*5);
end
xlabel('15分钟时间间隔');
ylabel('负荷数据');
legend(['原始数据', '场景数据 (概率:', num2str(probabilities(1)), ')'], ...
    ['场景数据 (概率:', num2str(probabilities(2)), ')'], ...
    ['场景数据 (概率:', num2str(probabilities(3)), ')'], ...
    ['场景数据 (概率:', num2str(probabilities(4)), ')'], ...
    ['场景数据 (概率:', num2str(probabilities(5)), ')']);

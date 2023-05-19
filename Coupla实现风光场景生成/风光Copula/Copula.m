function [C_solar,C_wind] = Copula(solardata,winddata)
    %��������
    solardata = solardata.*10;
    winddata = winddata.*10; 
    %�����ع�������ȡǰ92�������
    group_size = 96;
    each_season_days = 92;
    season_days = group_size*each_season_days;
    num_groups = season_days / group_size;
    
    winddata = reshape(winddata(1:season_days, 1),group_size,num_groups).';
    solardata = reshape(solardata(1:season_days, 1),group_size,num_groups).';
    
    scenarionum = 500;  % ��ʼ������Ŀ�����޸�
    num_cluster = 5;     % Ҫ�������ĳ�����Ŀ�����޸�
    ntime = 96;  % 96����
    
    % X��Y�ֱ�洢��͹��96��ʱ����ʷ�۲�����
    X = []; Y = [];
    for t = 1 : ntime
        X{t} = winddata(:, t);
        Y{t} = solardata(:, t);
    end
    %% Copula���
    % Frank-Copula ��������ͬʱ���Ǳ����ķǸ��븺��صĹ�ϵ
    % �ʲ��� Frank-Copula �����ֱ��96��ʱ�̽������
    
    for i = 1 : ntime
        U = ksdensity(X{i}, 'function', 'cdf'); % ���ܶȹ���
        V = ksdensity(Y{i}, 'function', 'cdf');
        
        alpha = copulafit('frank', [U(:) V(:)]); % ��ϳ��Ĳ���
        copulaparams.alpha = alpha;
        copulaparams.numParams = 1;
        copModels(i) = copulaparams;       
    end
    %% ���ƶ�ԪFrank-Copula���ܶȺ����ͷֲ�����ͼ
    [Udata, Vdata] = meshgrid(linspace(0,1,31));  % Ϊ��ͼ��Ҫ�������µ���������
    Ccdf_Frank = copulacdf('Frank', [Udata(:), Vdata(:)], copModels(48).alpha);
    
    
    %% ����
    Data = cell(1, ntime); 
    
    for j = 1: ntime
        Data{j} = copularnd('Frank', copModels(j).alpha, scenarionum); % ��Frank Copula�в���
    end
    
    %����任��ת��Ϊʵ�ʳ���
    w = cell(1, ntime);
    for k = 1 : ntime
        [fwind, xwind] = ecdf(X{k});
        funwind = @(pdwind)spline(fwind, xwind, pdwind);  % ������ֵ
        
        [fsolar, xsolar] = ecdf(Y{k});
        funsolar = @(pdsolar)spline(fsolar, xsolar, pdsolar);
        
        for j = 1 : scenarionum
            u = Data{k}(j, 1);
            v = Data{k}(j, 2);
            w{k}(j, 1) = max(funwind(u), 0);    % ȷ������ֵ�Ǹ�
            w{k}(j, 2) = max(funsolar(v), 0);
        end
    end
    %% ��������
    wind = [];
    solar = [];
    
    for i = 1 : ntime
        wind = [wind, w{i}(:, 1)];
        solar = [solar, w{i}(:, 2)];
    end
    
    % K-means����
    opts = statset('Display', 'final');
    [idx, C] = kmeans([wind, solar], num_cluster, 'Distance', 'cityblock', 'Replicates', 20, 'Options',opts);
    
    C_wind = C(:, 1:96);     % ���յķ�糡��
    C_solar = C(:, 97:end);  % ���յĹ������
    %% ���ʼ���
    for i = 1 : num_cluster
        p(i) =  length(find(idx == i)) / length(idx);
    end
    disp(['����������Ϊ: ', num2str(p)]);
    
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
    % ԭʼ��������ͼ
    figure(1);
    [ss,gg]=meshgrid(1:scenarionum,1:96);
    plot3(ss,gg,wind, 'linewidth', 1);
    title(['������������ɵķ�����', num2str(scenarionum), '������'])
    xlabel('����'); ylabel('ʱ��');zlabel('������ֵ');
    figure(2);
    plot3(ss,gg,solar,'linewidth', 1);
    title(['������������ɵĹ������', num2str(scenarionum), '������'])
    xlabel('����'); ylabel('ʱ��');zlabel('�������ֵ');
    % ��ͼ
    figure(3);
    [ss1,gg1]=meshgrid(1:num_cluster,1:96);
    plot3(ss1,gg1,C_wind, 'linewidth', 1);
    title(['�����������', num2str(num_cluster), '������'])
    xlabel('����'); ylabel('ʱ��');zlabel('������ֵ');
    figure(4);
    plot3(ss1,gg1,C_solar, 'linewidth', 1);
    title(['������������', num2str(num_cluster), '������'])
    xlabel('����'); ylabel('ʱ��');zlabel('�������ֵ');
    figure(5); 
    bar(p); 
    xlabel('������'); 
    ylabel('����');
    title('�������¸���');
    figure(6)
    plot(P_WD,'LineWidth',1.5); 
    xlabel('ʱ��/t'); 
    ylabel('����/kW');
    title('	��粻ȷ���Գ���');
    grid on
    figure(7)
    plot(P_PV,'LineWidth',1.5); 
    xlabel('ʱ��/t'); 
    ylabel('����/kW');
    title('	��粻ȷ���Գ���');
    grid on
    figure(8);  
    surf(Udata, Vdata, reshape(Ccdf_Frank,size(Udata)));  % ���ƶ�ԪFrank-Copula�ֲ�����ͼ
    xlabel('������'); 
    ylabel('�������'); 
    zlabel('C(u,v)');  
    title('��ԪFrank-Copula�ֲ�����ͼ')

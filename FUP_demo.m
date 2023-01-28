% Fast Unsupervised Dimensionality Reduction
% 快速无监督降维
% 
% min_{P>=0, P*1=1, W'StW=I, m}  sum（i=1，n）sum（j=1，m）pij|| W'xi-W'mj||^2 + r*||P||^2 
% 
% X  ：  d*n 原始样本空间, 每一列为一个样本
% P  ：  n*m 需要学习的到的 n 个样本和 m 个代表点之间的相似度矩阵
% W  ：  d*d1投影矩阵将样本空间从 d 维投影到 d1 维
% M  ：  d*m 代表点矩阵，从样本空间选取出来的 m 个代表点，用来构造 P ，实现“快速”的构图
% r  ：  正则化系数
% St ：  d*d 全局散度矩阵
%                      St = XHX'
%        其中，H=I-(1/n)*ones(n)
%        约束 W'StW=I 可以保证子空间数据在概率学上不相关
% Y  ：  d1*n 降维后的样本空间
% U  ：  d1*m 降维后的代表点矩阵
% m  ：  初始代表点的个数
% d1 ：  子空间维度
% K  ：  近邻数



clc
close all;
clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   需要修改数据集的名称     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Name = 'dig1-10_uni' ;          % 数据集名称
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([Name ,'.mat']);     % 加载数据集和标签
X = double(X);
X = full(X);
[n,~] = size(X);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      需要修改的参数         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 720;                % gamma 初始值
d1 = 9;                   % 子空间维度
c = 10;                   % 类别数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t = 100;                 % k-means 分类次数
lamda = 100;
iteration = 30;
K = 15;                   % FUP 近邻个数，在计算 gamma 时用到了
X = X(1:n,:)';
label = Y(1:n);
m = floor(0.8*n);         % FUP 初始代表点个数, 0.5 倍的样本个数


%% 初始化 M P
% 初始化 M ，采用 等距抽样 的策略
a = n/m;
idxx = a:a:n;
len = size(idxx,2);
if (len ~= m)
    error('代表点个数出错')
end
Idx = round(idxx);
M = X(:,Idx);

% 初始化 P ，自调整高斯方法 
distXM = Eu2_distance(X,M);                       % 计算每一个(样本)和(代表点)之间的欧式距离
distXM_sort = sort(distXM,2);                     % 每一行按升序排序
P = zeros(n,m);                                   % 创建 n*m 的矩阵 P
for i = 1:n
    for j = 1:m
        diK = (distXM_sort(i,8))^0.5;             % 第 i 个样本点距离它的最近的第 7 个近邻的距离
        djK = (distXM_sort(Idx(j),8))^0.5;        % 第 j 个代表点距离它的最近的第 7 个近邻的距离
        sigma = diK * djK + eps;
        P(i,j) = exp(-distXM(i,j)/sigma);   % n*m 初始化的 P ，采用自调整高斯方法构建
    end
end

objva = zeros(1,iteration);
H = diag(ones(1,n)) - (1/n)*ones(n);        % 定心矩阵
St =  X*H*X';                               % 全局散度矩阵

%% 迭代求解
for iter = 1:iteration
    % 固定 P M 求 W
    D_column = diag(sum(P,1));                  % 矩阵 P 的列和
    A = X*X'-(X*P*M'+M*P'*X')+M*D_column*M';    % 由于 MATLAB 计算的原因，这里计算出的 A 可能不是对称矩阵
    A = (A + A')/2;                             % 这一步确保A为对称矩阵
    St = (St + St')/2;
    [Geigve,Geigva] = eig(A,St);                % 求解 A 和 St 的广义特征值和广义特征向量
    Geigva = diag(Geigva);                      % 将特征值矩阵转换为特征值向量
    [~, index] = sort(Geigva);                  % 特征值升序排列
    eigve_sort = Geigve(:,index);               % 特征向量按照特征值的顺序对应排列
    W = eigve_sort(:,1:d1);                     % 选取前 d1 个最小的特征值对应的特征向量作为 W 
    
    % 固定 P W 求 M
    Y = W'*X;                                   % 降维后的样本空间
    for j = 1:m
        U(:,j) = (Y*P(:,j))/(D_column(j,j));    % 更新投影之后的代表点矩阵
    end
    distYU1 = Eu2_distance(Y,U);                % 计算投影后样本点和代表点之间的欧式距离
    [~,idx] = sort(distYU1,1);                  % 将样本点和代表点之间的距离进行排序
    
    for k = 1:size(idx,2)
        idx_near = idx(1:15,k);
        M_update(k,:) = mean(X(:,idx_near)');   % 选择距离代表点最近的 K 个样本点作为代表点的近邻，取平均作为 M
    end
    M = M_update';
    U = W'*M;
    m = size(M,2);                              % 更新代表点个数
    
    % 固定 W M 求 P
    distYU = Eu2_distance(Y,U);                 % 更新子空间内样本点和代表点之间的欧式距离矩阵
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     需要修改      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     distYU_sort = sort(distYU,2);               % 每一行按升序排序
%     r = (K*sum(distYU_sort(:,K+1)) - sum(sum(distYU_sort(:,1:K),1),2))/2/n + eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for i=1:n
        ad = -distYU(i,:)/2/r;
        P(i,:) = EProjSimplex_new(ad);          % 更新 P 
    end
    objva(1,iter) = trace(Y*Y'+U*diag(sum(P,1))*U'-2*Y*P*U'); % 每一步迭代之后目标函数值
end




%% 计算分类结果
[ Y2 ] = PCA(X,d1);               % PCA
[ Y3 ] = LPP(X,d1);               % LPP
Y1 = W'*X;                        % FUP
Y1 = real(Y1);
for q=1:t
    [idx1] = kmeans(Y1',c);       % 使用 K-Means 聚类
    [idx2] = kmeans(Y2',c);
    [idx3] = kmeans(Y3',c);
    R(q,1:3) = ClusteringMeasure(label,idx1);
    R(q,4:6) = ClusteringMeasure(label,idx2);
    R(q,7:9) = ClusteringMeasure(label,idx3);
end
FUP_result = zeros(2,3);
FUP_result(1,1) = mean(R(:,1));
FUP_result(2,1) = mean(R(:,2));
FUP_result(1,2) = mean(R(:,4));
FUP_result(2,2) = mean(R(:,5));
FUP_result(1,3) = mean(R(:,7));
FUP_result(2,3) = mean(R(:,8));
FUP_result
figure(1)
plot(1:iter,objva,'-s'),hold on
xlabel('Iteration times')
ylabel('Objective function value')
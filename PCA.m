function [ Y ] = PCA( X,d1 )
% 实现主成分分析（Principal Component Analysis， PCA）的函数文件
% 函数返回降维后的样本数据
% d1  ： 降维后的维度 
% X   ： d*n 每一列为一个样本
% W   ： d*d1 投影矩阵
% Y   ： d1*n 降维后的样本，每一列为一个样本

%% 数据中心化
X = X - mean(X')';

%% 对 X 进行奇异值分解
% [U,S,~] = svd(X);      % 对 X 进行奇异值分解，左奇异矩阵 U 为 XX' 的特征向量矩阵
[U,S] = eig(X*X');
S = diag(S);
[~,ind] = sort(S,'descend');  % 对奇异值进行 降序 排序
U_sort = U(:,ind);            % 按照奇异值的顺序对左奇异矩阵的列向量进行排列
W = U_sort(:,1:d1);           % 选取前 d1 个最大的奇异值对应的左奇异矩阵的列向量为 W
Y = W'*X;
end
function [ Y ] = LPP( X,d1 )
% 局部保持投影（Local Preserving Projection, LPP）
% min  sum_(i,j=1:n)[S_ij*(y_i-y_j)^2]
% 
% X  ： d*n 样本矩阵，每一列为一个样本
% Y  ： d1*n 投影后的样本矩阵
% W  ： d*d1 投影矩阵
% d1 ： 投影后的维度
% S  ： n*n 相似度矩阵，由自调整高斯方法构建
% K  ： 自调整高斯方法设定的近邻数

% 构建 S ，自调整高斯方法 
[~,n] = size(X);
distXX = Eu2_distance(X,X);                       % 计算每一个(样本)和(代表点)之间的欧式距离
distXX_sort = sort(distXX,2);                     % 每一行按升序排序
S = zeros(n,n);                                   % 创建 n*m 的矩阵 P
for i = 1:n
    for j = 1:n
        diK = (distXX_sort(i,8))^0.5;             % 第 i 个样本点距离它的最近的第 7 个近邻的距离
        djK = (distXX_sort(j,8))^0.5;             % 第 j 个代表点距离它的最近的第 7 个近邻的距离
        sigma = diK * djK + eps;
        S(i,j) = exp(-distXX(i,j)/sigma);         % n*n 采用自调整高斯方法构建的邻接矩阵 S
    end
end

D=diag(sum(S,2));                                 % 度矩阵
L=D-S;                                            % 拉普拉斯矩阵
A = X*L*X';
B = X*D*X';
[Gve,Gva] = eig(A,B);                             % 求解广义特征值特征向量问题
[~,ind] = sort(diag(Gva));
Gve_sort = Gve(:,ind);
W = Gve_sort(:,1:d1);                             % 取前 d1 个最小的特征值对应的特征向量为投影矩阵
if ~(isreal(W))
    W = abs(W);
end
Y = W'*X;                                         % 投影后的样本
end
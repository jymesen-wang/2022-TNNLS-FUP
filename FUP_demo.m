% Fast Unsupervised Dimensionality Reduction
% �����޼ල��ά
% 
% min_{P>=0, P*1=1, W'StW=I, m}  sum��i=1��n��sum��j=1��m��pij|| W'xi-W'mj||^2 + r*||P||^2 
% 
% X  ��  d*n ԭʼ�����ռ�, ÿһ��Ϊһ������
% P  ��  n*m ��Ҫѧϰ�ĵ��� n �������� m �������֮������ƶȾ���
% W  ��  d*d1ͶӰ���������ռ�� d άͶӰ�� d1 ά
% M  ��  d*m �������󣬴������ռ�ѡȡ������ m ������㣬�������� P ��ʵ�֡����١��Ĺ�ͼ
% r  ��  ����ϵ��
% St ��  d*d ȫ��ɢ�Ⱦ���
%                      St = XHX'
%        ���У�H=I-(1/n)*ones(n)
%        Լ�� W'StW=I ���Ա�֤�ӿռ������ڸ���ѧ�ϲ����
% Y  ��  d1*n ��ά��������ռ�
% U  ��  d1*m ��ά��Ĵ�������
% m  ��  ��ʼ�����ĸ���
% d1 ��  �ӿռ�ά��
% K  ��  ������



clc
close all;
clear all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ��Ҫ�޸����ݼ�������     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Name = 'dig1-10_uni' ;          % ���ݼ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([Name ,'.mat']);     % �������ݼ��ͱ�ǩ
X = double(X);
X = full(X);
[n,~] = size(X);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      ��Ҫ�޸ĵĲ���         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 720;                % gamma ��ʼֵ
d1 = 9;                   % �ӿռ�ά��
c = 10;                   % �����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t = 100;                 % k-means �������
lamda = 100;
iteration = 30;
K = 15;                   % FUP ���ڸ������ڼ��� gamma ʱ�õ���
X = X(1:n,:)';
label = Y(1:n);
m = floor(0.8*n);         % FUP ��ʼ��������, 0.5 ������������


%% ��ʼ�� M P
% ��ʼ�� M ������ �Ⱦ���� �Ĳ���
a = n/m;
idxx = a:a:n;
len = size(idxx,2);
if (len ~= m)
    error('������������')
end
Idx = round(idxx);
M = X(:,Idx);

% ��ʼ�� P ���Ե�����˹���� 
distXM = Eu2_distance(X,M);                       % ����ÿһ��(����)��(�����)֮���ŷʽ����
distXM_sort = sort(distXM,2);                     % ÿһ�а���������
P = zeros(n,m);                                   % ���� n*m �ľ��� P
for i = 1:n
    for j = 1:m
        diK = (distXM_sort(i,8))^0.5;             % �� i �������������������ĵ� 7 �����ڵľ���
        djK = (distXM_sort(Idx(j),8))^0.5;        % �� j ������������������ĵ� 7 �����ڵľ���
        sigma = diK * djK + eps;
        P(i,j) = exp(-distXM(i,j)/sigma);   % n*m ��ʼ���� P �������Ե�����˹��������
    end
end

objva = zeros(1,iteration);
H = diag(ones(1,n)) - (1/n)*ones(n);        % ���ľ���
St =  X*H*X';                               % ȫ��ɢ�Ⱦ���

%% �������
for iter = 1:iteration
    % �̶� P M �� W
    D_column = diag(sum(P,1));                  % ���� P ���к�
    A = X*X'-(X*P*M'+M*P'*X')+M*D_column*M';    % ���� MATLAB �����ԭ������������ A ���ܲ��ǶԳƾ���
    A = (A + A')/2;                             % ��һ��ȷ��AΪ�Գƾ���
    St = (St + St')/2;
    [Geigve,Geigva] = eig(A,St);                % ��� A �� St �Ĺ�������ֵ�͹�����������
    Geigva = diag(Geigva);                      % ������ֵ����ת��Ϊ����ֵ����
    [~, index] = sort(Geigva);                  % ����ֵ��������
    eigve_sort = Geigve(:,index);               % ����������������ֵ��˳���Ӧ����
    W = eigve_sort(:,1:d1);                     % ѡȡǰ d1 ����С������ֵ��Ӧ������������Ϊ W 
    
    % �̶� P W �� M
    Y = W'*X;                                   % ��ά��������ռ�
    for j = 1:m
        U(:,j) = (Y*P(:,j))/(D_column(j,j));    % ����ͶӰ֮��Ĵ�������
    end
    distYU1 = Eu2_distance(Y,U);                % ����ͶӰ��������ʹ����֮���ŷʽ����
    [~,idx] = sort(distYU1,1);                  % ��������ʹ����֮��ľ����������
    
    for k = 1:size(idx,2)
        idx_near = idx(1:15,k);
        M_update(k,:) = mean(X(:,idx_near)');   % ѡ�������������� K ����������Ϊ�����Ľ��ڣ�ȡƽ����Ϊ M
    end
    M = M_update';
    U = W'*M;
    m = size(M,2);                              % ���´�������
    
    % �̶� W M �� P
    distYU = Eu2_distance(Y,U);                 % �����ӿռ���������ʹ����֮���ŷʽ�������
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     ��Ҫ�޸�      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     distYU_sort = sort(distYU,2);               % ÿһ�а���������
%     r = (K*sum(distYU_sort(:,K+1)) - sum(sum(distYU_sort(:,1:K),1),2))/2/n + eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for i=1:n
        ad = -distYU(i,:)/2/r;
        P(i,:) = EProjSimplex_new(ad);          % ���� P 
    end
    objva(1,iter) = trace(Y*Y'+U*diag(sum(P,1))*U'-2*Y*P*U'); % ÿһ������֮��Ŀ�꺯��ֵ
end




%% ���������
[ Y2 ] = PCA(X,d1);               % PCA
[ Y3 ] = LPP(X,d1);               % LPP
Y1 = W'*X;                        % FUP
Y1 = real(Y1);
for q=1:t
    [idx1] = kmeans(Y1',c);       % ʹ�� K-Means ����
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
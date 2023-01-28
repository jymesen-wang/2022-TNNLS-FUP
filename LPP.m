function [ Y ] = LPP( X,d1 )
% �ֲ�����ͶӰ��Local Preserving Projection, LPP��
% min  sum_(i,j=1:n)[S_ij*(y_i-y_j)^2]
% 
% X  �� d*n ��������ÿһ��Ϊһ������
% Y  �� d1*n ͶӰ�����������
% W  �� d*d1 ͶӰ����
% d1 �� ͶӰ���ά��
% S  �� n*n ���ƶȾ������Ե�����˹��������
% K  �� �Ե�����˹�����趨�Ľ�����

% ���� S ���Ե�����˹���� 
[~,n] = size(X);
distXX = Eu2_distance(X,X);                       % ����ÿһ��(����)��(�����)֮���ŷʽ����
distXX_sort = sort(distXX,2);                     % ÿһ�а���������
S = zeros(n,n);                                   % ���� n*m �ľ��� P
for i = 1:n
    for j = 1:n
        diK = (distXX_sort(i,8))^0.5;             % �� i �������������������ĵ� 7 �����ڵľ���
        djK = (distXX_sort(j,8))^0.5;             % �� j ������������������ĵ� 7 �����ڵľ���
        sigma = diK * djK + eps;
        S(i,j) = exp(-distXX(i,j)/sigma);         % n*n �����Ե�����˹�����������ڽӾ��� S
    end
end

D=diag(sum(S,2));                                 % �Ⱦ���
L=D-S;                                            % ������˹����
A = X*L*X';
B = X*D*X';
[Gve,Gva] = eig(A,B);                             % ����������ֵ������������
[~,ind] = sort(diag(Gva));
Gve_sort = Gve(:,ind);
W = Gve_sort(:,1:d1);                             % ȡǰ d1 ����С������ֵ��Ӧ����������ΪͶӰ����
if ~(isreal(W))
    W = abs(W);
end
Y = W'*X;                                         % ͶӰ�������
end
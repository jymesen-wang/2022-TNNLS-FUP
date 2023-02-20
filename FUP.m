function [ W,M,P,objva,r ] = FUP( X,P,m,d1,K,iteration,a)

[~,n] = size(X);

if a == 1
    a = n/m;
    idxx = a:a:n;
    Idx = round(idxx);
    M = X(:,Idx);
else
    [~,M] = kmeans(X',m);
    M = M';
end

objva = zeros(1,iteration);
H = diag(ones(1,n)) - (1/n)*ones(n);    
St =  X*H*X';                             


for iter = 1:iteration
    D_column = diag(sum(P,1));                  
    A = X*X'-(X*P*M'+M*P'*X')+M*D_column*M';    
    A = (A + A')/2;                            
    St = (St + St')/2;
    [Geigve,Geigva] = eig(A,St);               
    Geigva = diag(Geigva);                    
    [~, index] = sort(Geigva);                 
    eigve_sort = Geigve(:,index);               
    W = eigve_sort(:,1:d1);                     
    

    Y = W'*X;                                 
    M = X*P*pinv(D_column);
    U = W'*M;
    

    distYU = Eu2_distance(Y,U);                 
    distYU_sort = sort(distYU,2);             
    r = (K*sum(distYU_sort(:,K+1)) - sum(sum(distYU_sort(:,1:K),1),2))/2/n + eps;
    for i=1:n
        ad = -distYU(i,:)/2/r;
        P(i,:) = EProjSimplex_new(ad);        
    end
    objva(1,iter) = trace(Y*Y'+U*diag(sum(P,1))*U'-2*Y*P*U'); 
end
end
% 子空间为2维，进行类别可视化
function [] = ClassVisualization(idx1,idx2,idx3,Y1,Y2,Y3,n)
    figure(1)
    for j = 1:n

%                 subplot(1,4,1); pbaspect([1 1 1]);  % 原始数据的结果
%                 if (label(j)==1)
%                     plot(X(1,j),X(2,j),'bo'),hold on
%                 elseif (label(j)==2)
%                     plot(X(1,j),X(2,j),'r^'),hold on
%                 elseif (label(j)==3)
%                     plot(X(1,j),X(2,j),'g*'),hold on
%                 end

        subplot(1,3,1); pbaspect([1 1 1]);  % FUP的结果
        if (idx1(j)==1)
            plot(Y1(1,j),Y1(2,j),'bo','MarkerSize',10),hold on
        elseif (idx1(j)==2)
            plot(Y1(1,j),Y1(2,j),'r^','MarkerSize',10),hold on
        elseif (idx1(j)==3)
            plot(Y1(1,j),Y1(2,j),'g*','MarkerSize',10),hold on
        elseif (idx1(j)==4)
            plot(Y1(1,j),Y1(2,j),'mx','MarkerSize',10),hold on
        elseif (idx1(j)==5)
            plot(Y1(1,j),Y1(2,j),'ks','MarkerSize',10),hold on
        end

        subplot(1,3,2);  pbaspect([1 1 1]); % PCA 的结果
        if (idx2(j)==1)
            plot(Y2(1,j),Y2(2,j),'bo','MarkerSize',10),hold on
        elseif (idx2(j)==2)
            plot(Y2(1,j),Y2(2,j),'r^','MarkerSize',10),hold on
        elseif (idx2(j)==3)
            plot(Y2(1,j),Y2(2,j),'g*','MarkerSize',10),hold on
        elseif (idx2(j)==4)
            plot(Y2(1,j),Y2(2,j),'mx','MarkerSize',10),hold on
        elseif (idx2(j)==5)
            plot(Y2(1,j),Y2(2,j),'ks','MarkerSize',10),hold on
        end

        subplot(1,3,3);  pbaspect([1 1 1]); % LPP 的结果
        if (idx3(j)==1)
            plot(Y3(1,j),Y3(2,j),'bo','MarkerSize',10),hold on
        elseif (idx3(j)==2)
            plot(Y3(1,j),Y3(2,j),'r^','MarkerSize',10),hold on
        elseif (idx3(j)==3)
            plot(Y3(1,j),Y3(2,j),'g*','MarkerSize',10),hold on
        elseif (idx3(j)==4)
            plot(Y3(1,j),Y3(2,j),'mx','MarkerSize',10),hold on
        elseif (idx3(j)==5)
            plot(Y3(1,j),Y3(2,j),'ks','MarkerSize',10),hold on
        end
    end
%             J(1) = subplot(1,3,1);
    J(2) = subplot(1,3,1);
    J(3) = subplot(1,3,2);
    J(4) = subplot(1,3,3);

%             title( J(1),'Original data','FontSize',15)
    title( J(2),'FUP','FontSize',15)
    title( J(3),'PCA','FontSize',15)
    title( J(4),'LPP','FontSize',15)
end
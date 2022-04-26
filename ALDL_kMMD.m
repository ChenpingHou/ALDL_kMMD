function [W,a,obj] = ALDL_kMMD(X1,X2,y,option)
%% 参数设定
[d,n1]=size(X1);%d为数据维数，n1为标签数据数
[d,n2]=size(X2);%n2为未标签数据数
[c,n1] = size(y);%c为标签类数
t = option.t;%第二项正则化项的系数
b = option.b;%选点个数
iter  = option.iter;
a0 = b/n2*ones(1,n2);
a = a0;
for v = 1:iter
    %% 更新W
    %    T1 = zeros();
    T1 = (1/(n1+b))*sum(X1,2);
    T2 = 0;
    for i = 1:n2
        T2 = T2 + a(i)*X2(:,i);
    end
    T3 = (1/(n1+b))*T2;
    T4 = 0;
    T4 = (1/(n2-b))*(sum(X2,2)-T2);
    T = T1+T3-T4;
    W = y*X1'/(X1*X1'+ t*T*T');
    %% 更新a
    X = [X2,X1];
    %     K = zeros(n1+n2,n1+n2);
    K = (W*X)'*(W*X);
    K1 = zeros(n2,n2);
    for i = 1:n2
        for j =1:n2
            K1(i,j) = K(i,j);
        end
    end
    K2 = zeros(n2,1);
    K22 = 0;
    for i =1:n2
        for j = 1:n2
            K22 = K22 +K1(i,j);
        end
        K2(i) = (n1+b)/(n1+n2)* K22;
        K22 = 0;
    end
    K3 = zeros(n2,1);
    K33 = 0;
    for i =1:n2
        for j = 1:n1
            K33 = K33 + K(i,n2+j);
        end
        K3(i) = (n2-b)/(n1+n2)* K33;
        K33 = 0;
    end
    K4 = -K2+K3;
    Aeq = ones(1,n2);
    beq = b;
    lb = zeros(1,n2);
    ub = ones(1,n2);
    a = quadprog(K1,K4,[],[],Aeq,beq,lb,ub);
    %    a0 = a;
    %% 计算obj
    T1 = zeros();
    T1 = (1/(n1+b))*sum(X1,2);
    T2 = 0;
    for i = 1:n2
        T2 = T2 + a(i)*X2(:,i);
    end
    T3 = (1/(n1+b))*T2;
    T4 = 0;
    T4 = (1/(n2-b))*(sum(X2,2)-T2);
    T = T1+T3-T4;
    
    obj(v) = (norm(W*X1-y,'fro'))^2 + t*(norm(W*T,'fro'))^2;
    if v>1
        if (obj(v-1) - obj(v))/obj(v)<1e-8
            break
        end
    end
end
    

function nearestpoly 
% compute the nearest polynomial to multiple given polynomials
% f0=1
% polynomial basis: e_j=x^j-1, j=1,2,...,n
% multiple given polynoial: f_i=x^i, i=1,2,...,m
% given zero: z=1;
% p=2,q, 0<q<=2
% written by Wenyu Hu
% date: Jul. 2019

clear all
clc
format long

q=.25;
n=6;
m=5;

z=2;

epsilon=1e-8;
isConverged=0;
% e_j(x)
Ve=zeros(n,1);
for j=1:n
    Ve(j)=z^(j-1);
end
norme=norm(Ve,2);

%
f0=@(x)1;
Vf0=zeros(m,1);
for i=1:m
    Vf0(i)=f0(z);
end

%
Vf=zeros(m,1);
for i=1:m
    f{i}=@(x)x^i;
    Vf(i)=f{i}(z);
end

%
if m<=n-1
    % initialize
    A=zeros(m,n);
    for i=1:m
        A(i,1)=-1;
        A(i,i+1)=1;
    end
    D=eye(m);
    X=zeros(m,n);
    c=zeros(n,1);
    iter=0;
tic    
    while ~isConverged
        iter=iter+1;
        % update c
        temp=(ones(m,1)'*D*Vf)/(norme^2);
        c=A'*D*ones(m,1)-temp*Ve;
        c=c/(ones(m,1)'*D*ones(m,1));
        
        % update X 
        pre_X=X;
        X=ones(m,1)*ones(m,1)'*D-(ones(m,1)'*D*ones(m,1))*eye(m);
        X=X*A-temp*ones(m,1)*Ve';
        X=X/(ones(m,1)'*D*ones(m,1));
        
        % update D
        epsilon1=1e-5;
        D=zeros(m);
        for i=1:m
            if norm(X(i,:),2)==0
                temp=norm(X(i,:),2)+epsilon1;
            else
                temp=norm(X(i,:),2);
            end  
            temp=temp^(2-q);
            D(i,i)=1/temp;
        end
        
        % Decision
        if iter>=1
            normX(iter)=Mixed_norm(X,q);
            gapX(iter)=abs(Mixed_norm(pre_X,q)-Mixed_norm(X,q))/Mixed_norm(X,q);
            if gapX(iter)<=epsilon
                isConverged=1;
            end
        end
    end
    
    TimeA=toc;
    
    disp(['Output: q=' num2str(q)])
    c
    % X
    iter
    normX
    gapX
    TimeA
    
    % f0(z)+c'*Ve
    
    figure(1); plot(1:iter,normX,'r-d','LineWidth',1.5); xlabel('iteration number');ylabel('||X_k||_{2,q}^q');
    figure(2); plot(2:iter,gapX(2:end),'b-s','LineWidth',1.5);xlabel('iteration number');ylabel('\rho_k');
else 
    disp('Error: m>=n');
end



%%%%%%%%%%%%%
%%%%%%%%%%%%%
function s=Mixed_norm(X,p)
% compute \|X\|_2,p^p

[m,n]=size(X);

s=0;
for k=1:m
    s=s+norm(X(k,:),2)^p;
end


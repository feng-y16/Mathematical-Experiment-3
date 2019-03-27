function [result_J,result_G,t] = compare_Jacobi_GS(A,b,m,tol,x0,x_gt)
t=zeros(2,1);
tic
[~,k_J,X_J]=Jacobi_inf(A,b,m,tol,x0);
t(1,1)=toc;
tic
[~,k_G,X_G]=GS_inf(A,b,m,tol,x0);
t(2,1)=toc;
result_J=zeros(k_J+1,1);
for i=1:1:k_J+1
    result_J(i)=norm(X_J(:,i)-x_gt,inf);
end
result_G=zeros(k_G+1,1);
for i=1:1:k_G+1
    result_G(i)=norm(X_G(:,i)-x_gt,inf);
end
end


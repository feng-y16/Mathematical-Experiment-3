function [x,k,X,error]=Jacobi_output_error(A,b,tol,x0,x_gt)
D=diag(diag(A));  
n=length(A);
L=-(tril(A)-D);  
U=-(triu(A)-D);
f=D\b;   
B=D\(L+U);
x=x0*ones(n,1);
X=x;
k=1; 
x=B*x+f;
X=[X,x];
while norm(X(:,k+1)-X(:,k),inf)>tol
    k=k+1;
    x=B*x+f;
    X=[X,x];
end
error=zeros(k+1,1);
for i=1:1:k+1
    error(i)=norm(X(:,i)-x_gt,inf);
end
end

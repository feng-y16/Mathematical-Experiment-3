function [x,k,X]=Jacobi_inf(A,b,m,tol,x0,varargin)
D=diag(diag(A));  
n=length(A);
L=-(tril(A)-D);  
U=-(triu(A)-D);
f=D\b;   
B=D\(L+U);
if nargin==4
    x=zeros(n,1); 
else 
    x=x0*ones(n,1);
end
X=x;
k=0; 
while norm(A*x-b,inf)/norm(b,inf)>tol && k<m
    k=k+1;
    x=B*x+f;
    X=[X,x];
end

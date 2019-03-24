function [x,k,X]=GS(A,b,m,tol,x0,varargin)
n=length(A);  
D=diag(diag(A));
L=tril(A);  
U=-(triu(A)-D);
f=L\b;   
B=L\U;
if nargin==4
    x=zeros(n,1); 
else 
    x=x0*ones(n,1);
end
X=x;
k=0;
while norm(A*x-b)/norm(b)>tol && k<m
    k=k+1;
    x=B*x+f;
    X=[X,x];
end


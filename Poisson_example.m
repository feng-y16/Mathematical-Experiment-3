function [error,t]=Poisson_example(n,name,f,f_m,tol,varargin)
tic
m=n*n;  
c=linspace(0,1,n+2); 
h=1/(n+1);
A=4*eye(m)+diag(-ones(m-1,1),1)+diag(-ones(m-1,1),-1)+diag(-ones(m-n,1),n)+diag(-ones(m-n,1),-n);
for k=2:n-1
    A((k-1)*n+1,(k-1)*n)=0; 
    A(k*n,k*n+n)=0;
end
A(n*(n-1)+1,n*(n-1))=0;
A(n,n+n)=0;

xx=c(2:n+1)'*ones(1,n); 
yy=ones(n,1)*c(2:n+1); 
zz=2*pi.^2*sin(pi*xx).*sin(pi.*yy);
b=h^2*reshape(zz,m,1);
if nargin>2
    r=f(A,b,f_m,tol);
else
    r=A\b; 
end
r=reshape(r,n,n); 
r=[zeros(1,n+2);zeros(n,1),r,zeros(n,1);zeros(1,n+2)];
t=toc;
r_gt=zeros(n+2,n+2);
for i=1:1:n+2
    for j=1:1:n+2
        r_gt(i,j)=sin(pi*c(i))*sin(pi*c(j));
    end
end
error=sum(sum(abs(r_gt-r)))/(n+2)^2;

[X,Y]=meshgrid(c, c);  
mesh(X,Y,r,'LineWidth',3);
set(gca,'FontSize',28)
xlabel('x')
ylabel('y')
zlabel('u')
title('Poisson example')


set(gcf,'outerposition',get(0,'screensize'));
saveas(gcf,['3_' num2str(n) '_' name '.png'])
close
end

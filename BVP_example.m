function [error,t]=BVP_example(num,name,f,f_m,tol,varargin)
tic
h=1.0/num;
x=1:h:2;
y0=1/2; 
yn=4+4*log(2);
n=length(x); 
n=n-1;
p=2./x; 
q=-6./(x.^2);
r=5-6*x+7*x.^2;

b=-2+h^2*q(2:n); 
a=1-h*p(2:n)/2; 
c=1+h*p(2:n)/2;
d=h^2*r(2:n);
d(1)=d(1)-a(1)*y0; 
d(n-1)=d(n-1)-c(n-1)*yn;

A=diag(b)+diag(a(2:n-1),-1)+diag(c(1:n-2),1);

if nargin>2
    y_temp=f(A,d',f_m,tol);
else
    y_temp=A\d'; 
end

y=[y0; y_temp; yn];  
t=toc;
z=x.^2-x.^3+(x.^4)/2+(x.^2).*log(x);
error=sum(abs(z'-y))/n;

plot(x,y,'LineWidth',3);
hold
plot(x,z,'--','LineWidth',3)
set(gca,'FontSize',28)
xlabel('x')
ylabel('y')
legend('估计曲线','真实曲线')
title('BVP example')

set(gcf,'outerposition',get(0,'screensize'));
saveas(gcf,['3_2_' num2str(n) '_' name '.png'])
close

end
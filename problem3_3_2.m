tol=1.7e-3;
m=20000;
t=zeros(6,1);

[error11,t(1)]=Poisson_example(30,'');
[error12,t(2)]=Poisson_example(30,'Jacobi',@Jacobi,m,tol);
[error13,t(3)]=Poisson_example(30,'GS',@GS,m,tol);

tol=1e-8;
[error21,t(4)]=Poisson_example(60,'');
[error22,t(5)]=Poisson_example(60,'Jacobi',@Jacobi,m,tol);
tol=4.3e-4;
[error23,t(6)]=Poisson_example(60,'GS',@GS,m,tol);
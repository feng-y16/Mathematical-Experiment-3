tol=1.7e-3;
m=20000;
t=zeros(6,1);

tic
error11=Poisson_example(30,'');
t(1)=toc;
tic
error12=Poisson_example(30,'Jacobi',@Jacobi,m,tol);
t(2)=toc;
tic
error13=Poisson_example(30,'GS',@GS,m,tol);
t(3)=toc;

tol=1e-8;
tic
error21=Poisson_example(60,'');
t(4)=toc;
tic
error22=Poisson_example(60,'Jacobi',@Jacobi,m,tol);
t(5)=toc;
tol=4.3e-4;
tic
error23=Poisson_example(60,'GS',@GS,m,tol);
t(6)=toc;
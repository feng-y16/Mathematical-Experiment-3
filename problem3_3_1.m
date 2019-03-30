m=20000;
t=zeros(6,1);

[error11,t(1)]=BVP_example(100,'');
tol=1.78e-7;
[error12,t(2)]=BVP_example(100,'Jacobi',@Jacobi,m,tol);
tol=1.02e-7;
[error13,t(3)]=BVP_example(100,'GS',@GS,m,tol);

[error21,t(4)]=BVP_example(1000,'');
tol=1e-11;
[error22,t(5)]=BVP_example(1000,'Jacobi',@Jacobi,m,tol);
tol=4.3e-11;
[error23,t(6)]=BVP_example(1000,'GS',@GS,m,tol);
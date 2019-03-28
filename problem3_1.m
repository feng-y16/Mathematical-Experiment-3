n=5;
t1=zeros(3,1);
[V,b_v] = Vandermonde_and_b(n);
[H,b_h] = Hilbert_and_b(n);
tic
error_v=norm(V\b_v-ones(n,1),2);
error_h=norm(H\b_h-ones(n,1),2);
t1(1)=toc;
tic
error_v_Jacobi=norm(Jacobi(V,b_v,1e+100,error_v,0)-ones(n,1),2);
error_h_Jacobi=norm(Jacobi(H,b_h,1e+100,error_h,0)-ones(n,1),2);
t1(2)=toc;
tic
error_v_GS=norm(GS(V,b_v,1e+100,error_v*1e-10,0)-ones(n,1),2);
error_h_GS=norm(GS(H,b_h,1e+100,error_h*1e-10,0)-ones(n,1),2);
t1(3)=toc;

epsilon=[1e-10,1e-8,1e-6];
tol=1e-5;
m=2000;
x0=0;

t2=zeros(4,1);
tic
data=collect_perturb_data(@purturb,epsilon,4);
t2(1)=toc;
tic
data_Jacobi=collect_perturb_data(@purturb_Jacobi,epsilon,4,tol,m,x0);
t2(2)=toc;
tic
data_GS=collect_perturb_data(@purturb_GS,epsilon,4,tol,m,x0);
t2(3)=toc;
tic
data_PCG=collect_perturb_data(@purturb_PCG,epsilon,4,tol,m,x0);
t2(4)=toc;

norm_Jacobi_GS=compute_norm(4);
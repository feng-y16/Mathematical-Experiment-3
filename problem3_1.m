n=5;
[V,b_v] = Vandermonde_and_b(n);
[H,b_h] = Hilbert_and_b(n);
error_v=norm(V\b_v-ones(n,1),2);
error_h=norm(H\b_h-ones(n,1),2);
fprintf("error_v=%f\n",error_v);
fprintf("error_h=%f\n",error_h);

epsilon=[1e-10,1e-8,1e-6];
tol=1e-5;
m=2000;
x0=0;

t=zeros(4,1);
tic
data=collect_perturb_data(@purturb,epsilon,4);
t(1)=toc;
tic
data_Jacobi=collect_perturb_data(@purturb_Jacobi,epsilon,4,tol,m,x0);
t(2)=toc;
tic
data_GS=collect_perturb_data(@purturb_GS,epsilon,4,tol,m,x0);
t(3)=toc;
tic
data_PCG=collect_perturb_data(@purturb_PCG,epsilon,4,tol,m,x0);
t(4)=toc;

norm_Jacobi_GS=compute_norm(4);
A_1=generate_A(1);
tol=1e-5;
m=2000;
t=zeros(4,1);

b_1=sum(transpose(A_1))';
x_gt1=ones(20,1);
x0=0;
tic
[J1,G1,t(1:2)]=compare_Jacobi_GS(A_1,b_1,m,tol,x0,x_gt1);
x0=-10;
[J2,G2,t(3:4)]=compare_Jacobi_GS(A_1,b_1,m,tol,x0,x_gt1);

x_gt2=random('Normal',5,5,20,1);
b_2=A_1*x_gt2;
x0=0.1;
[J3,G3,]=compare_Jacobi_GS(A_1,b_2,m,tol,x0,x_gt2);
x0=-10;
[J4,G4,]=compare_Jacobi_GS(A_1,b_2,m,tol,x0,x_gt2);

i_num=8;
k=zeros(i_num,1);
error=zeros(i_num,50);
for i=1:1:i_num
    multiplier=2^i;
    A=generate_A(multiplier);
    b=sum(transpose(A))';
    x_gt=ones(20,1);
    x0=0;
    [~,k(i),~,error_temp]=Jacobi_output_error(A,b,tol,x0,x_gt);
    error_temp=error_temp';
    len_error=length(error_temp);
    error(i,1:len_error)=error_temp;
end
    
    


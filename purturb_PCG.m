function [data] = purturb_PCG(epsilon,n_num,tol,m,x0,varargin)
data=zeros(n_num,10);
for i=1:1:n_num
    n=5+2*(i-1);
    [V,b_v] = Vandermonde_and_b(n);
    [H,b_h] = Hilbert_and_b(n);
    data(i,1)=cond(V,2);
    data(i,2)=cond(H,2);
    %fprintf("cond(V_%d)=%f\n",n,data(i,1))
    %fprintf("cond(H_%d)=%f\n",n,data(i,2))
    
    V(n,n)=V(n,n)+epsilon;
    H(n,n)=H(n,n)+epsilon;
    if nargin==4
        data(i,3)=norm(pcg(V,b_v,tol,m)-ones(n,1),2)/norm(ones(n,1),2);
        data(i,4)=norm(pcg(H,b_h,tol,m)-ones(n,1),2)/norm(ones(n,1),2);
    else
        data(i,3)=norm(pcg(V,b_v,tol,m,[],[],x0*ones(n,1))-ones(n,1),2)/norm(ones(n,1),2);
        data(i,4)=norm(pcg(H,b_h,tol,m,[],[],x0*ones(n,1))-ones(n,1),2)/norm(ones(n,1),2);
    end
    %fprintf("error_v=%f\n",data(i,3));
    %fprintf("error_h=%f\n",data(i,4));
    
    V(n,n)=V(n,n)-epsilon;
    H(n,n)=H(n,n)-epsilon;
    b_v(n)=b_v(n)+epsilon;
    b_h(n)=b_h(n)+epsilon;
    if nargin==4
        data(i,5)=norm(pcg(V,b_v,tol,m)-ones(n,1),2)/norm(ones(n,1),2);
        data(i,6)=norm(pcg(H,b_h,tol,m)-ones(n,1),2)/norm(ones(n,1),2);
    else
        data(i,5)=norm(pcg(V,b_v,tol,m,[],[],x0*ones(n,1))-ones(n,1),2)/norm(ones(n,1),2);
        data(i,6)=norm(pcg(H,b_h,tol,m,[],[],x0*ones(n,1))-ones(n,1),2)/norm(ones(n,1),2);
    end
    %fprintf("error_v=%f\n",data(i,5));
    %fprintf("error_h=%f\n",data(i,6));
    
    b_v(n)=b_v(n)-epsilon;
    b_h(n)=b_h(n)-epsilon;
    data(i,7)=data(i,1)*epsilon/norm(V,2);
    data(i,8)=data(i,2)*epsilon/norm(H,2);
    data(i,9)=data(i,1)*epsilon/norm(b_v,2);
    data(i,10)=data(i,2)*epsilon/norm(b_h,2);
end
end
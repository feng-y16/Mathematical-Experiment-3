function [data] = collect_perturb_data(f,epsilon,n_num,tol,m,x0,varargin)
n=length(epsilon);
data=zeros(4*n,10);
for i=1:1:n
    if nargin==3
        data(4*(i-1)+1:4*(i-1)+4,:)=f(epsilon(i),n_num);
    else
        if nargin==5
            data(4*(i-1)+1:4*(i-1)+4,:)=f(epsilon(i),n_num,tol,m);
        else
            data(4*(i-1)+1:4*(i-1)+4,:)=f(epsilon(i),n_num,tol,m,x0);
        end
    end
end
end


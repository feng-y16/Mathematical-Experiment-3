function [norm_Jacobi_GS] = compute_norm(n_num)
norm_Jacobi_GS=zeros(4,n_num);
for i=1:1:n_num
    n=5+2*(i-1);
    [V,] = Vandermonde_and_b(n);
    [H,] = Hilbert_and_b(n);
    D=diag(diag(V));
    L=-(tril(V)-D);  
    U=-(triu(V)-D);  
    norm_Jacobi_GS(1,i)=norm(D\(L+U),2);
    D=diag(diag(H));
    L=-(tril(H)-D);  
    U=-(triu(H)-D);  
    norm_Jacobi_GS(2,i)=norm(D\(L+U),2);
    D=diag(diag(V));
    L=tril(V);  
    U=-(triu(V)-D);   
    norm_Jacobi_GS(3,i)=norm(L\U,2);
    D=diag(diag(H));
    L=tril(H);  
    U=-(triu(H)-D);   
    norm_Jacobi_GS(4,i)=norm(L\U,2);
end


function [V,b] = Vandermonde_and_b(n)
V=zeros(n);
for i=1:1:n
    x=1+0.1*(i-1);
    for j=1:1:n
        V(i,j)=x^(j-1);
    end
end
b=sum(transpose(V))';
end
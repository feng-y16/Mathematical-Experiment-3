function [H,b] = Hilbert_and_b(n)
H=hilb(n);
b=sum(transpose(H))';
end
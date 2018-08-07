function s = FOCUSS(y,Phi)
iterations = 2;
N = size(Phi,2);
x = ones(N,1);

for it = 1:iterations
W_pk = diag(x);
matrix = ctranspose(Phi*W_pk);
A_plus = matrix / (Phi*W_pk * matrix);
q_k = A_plus*y;
x = W_pk*q_k;
end
s = abs(x).^2;
end
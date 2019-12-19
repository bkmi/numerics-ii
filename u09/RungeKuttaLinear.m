function [x, t, k] = RungeKuttaLinear(M, x0, I, tau, b, A)

% reformulation:
% k_i = M*g_i ,  i=1,...s
% z_i = g_i -x , i=1,...,s
% Z = [z_1;...;z_s].

b = b(:);
d = length(x0);
s = length(b);
t = I(1):tau:I(2);
N = length(t);
x = zeros(d, N);
k = zeros(d,s,N);
x(:,1) = x0;
DF = kron(A, M);
DG0 = eye(s*d) - tau*DF;

for i=2:N
    X = kron(ones(s,1), x(:,i-1));
    G0 = tau*DF*X;
    Z = DG0\G0;
    k(:,:,i) = M*reshape(X+Z, d, s);
    x(:,i) = x(:,i-1) + tau*k(:,:,i)*b;
end


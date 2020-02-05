a = [1, -0.1; 0.25, 1.1];
A = a' * a;
b = [1; 2];
u0 = [0.5; 0.4];
tol = 0.01;
uexact = A \ b;

[u, uk] = cg(A, b, u0, tol, uexact);
maxit = 1;
x = pcg(A,b,tol,maxit,eye(2),eye(2),u0);
printf(num2str(uk(:,1) == x))
printf("\n")


B = [1, 0.1; 0.1, 1];
Binv = inv(B);
pre = @(x)Binv * x;

cond = @(x)eig(x)(2)/eig(x)(1);
printf("condition numbers \n")
printf("kappa(A) = ")
printf(num2str(cond(A)))
printf("\n")
printf("kappa(B^{-1}A) = ")
printf(num2str(cond(pre(A))))
printf("\n")

[u_pre, uk_pre] = pcg_ben(A, b, u0, tol, uexact, pre);
x = pcg(A,b,tol,maxit,B,eye(2),u0);
printf(num2str(uk_pre(:,1) == x))
printf("\n")
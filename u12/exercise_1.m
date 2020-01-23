r0 = 0.01;
m = 0.1;
g = 1.62;
I = [0.0, 20.0];
p0 = 0.;
q0 = 0.;
taus = [1, 0.1];

for j=1:2
  figure;
  [p, q, t] = stoermer_verlet(m, g, r0, p0, q0, I, taus(j))
  q_mod = mod(q, 360.)
  plot(t, q_mod)
  title(strcat("Time evolution, Tau = ", num2str(taus(j))))
  xlabel("time")
  ylabel("q")
endfor

for j=1:2
  figure;
  plot(q, p)
  title(strcat("Phase Space evolution, Tau = ", num2str(taus(j))))
  xlabel("q")
  ylabel("p")
endfor

x = r0 * cos(q)
y = r0 * sin(q)

for j=1:2
  figure;
  plot(x, y)
  title(strcat("Euclidean Coordinates, Tau = ", num2str(taus(j))))
  xlabel("x")
  ylabel("y")
endfor

H = q + p

for j=1:2
  figure;
  plot(t, H)
  title(strcat("Hamiltonian in time, Tau = ", num2str(taus(j))))
  xlabel("t")
  ylabel("H")
endfor
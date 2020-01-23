r0 = 0.01;
m = 0.1;
g = 1.62;
I = [0.0, 20.0];
p0 = 0.;
q0 = 0.;
taus = [1, 0.1];
tau_names = {"1", "0p1"};

for j=1:2
  figure;
  [p, q, t] = stoermer_verlet(m, g, r0, p0, q0, I, taus(j));
  q_mod = mod(q, 360.);
  plot(t, q_mod);
  title(strcat("Time evolution, Tau = ", num2str(taus(j))));
  xlabel("time");
  ylabel("q");
  h = gcf;
  saveas(h, strcat("time_evo_tau_", tau_names{j}), "png")
endfor

for j=1:2
  figure;
  plot(q, p);
  title(strcat("Phase Space evolution, Tau = ", num2str(taus(j))));
  xlabel("q");
  ylabel("p");
  h = gcf;
  saveas(h, strcat("phase_space_evo_", tau_names{j}), "png")
endfor

x = r0 * cos(q);
y = r0 * sin(q);

for j=1:2
  figure;
  plot(x, y);
  title(strcat("Euclidean Coordinates, Tau = ", num2str(taus(j))));
  xlabel("x");
  ylabel("y");
  h = gcf;
  saveas(h, strcat("euc_evo_", tau_names{j}), "png")
endfor

H = q + p;

for j=1:2
  figure;
  plot(t, H);
  title(strcat("Hamiltonian in time, Tau = ", num2str(taus(j))));
  xlabel("t");
  ylabel("H");
  h = gcf;
  saveas(h, strcat("hamiltonian_evo_", tau_names{j}), "png")
endfor
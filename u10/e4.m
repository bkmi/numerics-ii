E = @(cj,c0) [cj, 0, -cj, 0, 0; 
              0, c0, 0, -c0, 0;
              -cj, 0, 2*cj, -cj, 0;
              0, -c0, -cj, c0+cj, 0;
              0, 0, 0, 0, 0];
A = @(G) -1 .* [G(1), 0, 0, 0, 0; 
          0, G(2)+G(4), 0, -G(4), 0;
          0, 0, G(3), 0, 0;
          0, -G(4), 0, G(4), 0;
          0, 0, 0, 0, G(5)]; 
b1 = @(alpha, g, u) [(1-alpha) * g(u(1) - u(3));
                     alpha * g(u(1) - u(3));
                     -g(u(1) - u(3)) - g(u(4) - u(3));
                     (1-alpha) * g(u(4) - u(3));
                     alpha * g(u(4) - u(3))];
b2 = @(G, Vin, Vdd)(t) -1 .* [G(1) * Vin(t);
                              G(2) * Vdd;
                              0;
                              0;
                              G(5) * Vdd];
b = b1 + b2;
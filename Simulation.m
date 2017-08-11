IC = [0 0 0 0 0 0 0 0 0 10];
opt = odeset('reltol', 1e-8, 'abstol', 1e-8);
[Time, X] = ode45(@DerivativesL, [0 5], IC, opt);


Animate(X);
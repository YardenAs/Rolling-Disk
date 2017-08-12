IC = [0 0 0 2.5*0.5 0 0 10*pi/180 0 0 2.5]; % Check that the initial condition satisfies the constraint
opt = odeset('reltol', 1e-8, 'abstol', 1e-8);
[Time, X] = ode45(@DerivativesL, [0 15], IC, opt);

Animate(Time, X);
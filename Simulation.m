IC = [0 0 0 4*0.5 0 0 10*pi/180 0 0 4]; 
% Check that the initial condition satisfies the no slip constraint
opt = odeset('reltol', 1e-8, 'abstol', 1e-8);
[Time, X] = ode45(@DerivativesL, [0 15], IC, opt);

Animate(X);
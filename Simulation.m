%% Lagrange and Newton derivation example comparison %%

IC = [0 0 0 10*0.5 0 0 10*pi/180 0 0 10]; 
% Check that the initial condition satisfies the no slip constraint
opt = odeset('reltol', 1e-8, 'abstol', 1e-8);
[TimeL, XL] = ode45(@DerivativesL, [0 15], IC, opt);
[TimeN, XN] = ode45(@DerivativesN, [0 15], IC, opt);
plot(TimeL, XL(:,7),TimeN, XN(:,7),'--r');
xlabel('Time [sec]'); ylabel('\theta [rad]');
legend('Lagrange Derivation', 'Newton Derivation');
% Animate(X);

%% Stability as function of dPhi %%
global theta0
global dtheta0
dtheta0 = 0;
theta0 = 15*pi/180;
dphi = linspace(theta0/5, theta0, 40)*180/pi;
% opt = odeset('reltol', 1e-9, 'abstol', 1e-9, 'Events', @isStable);
for ii = length(dphi):-1:1
    IC = [0 0 0 dphi(ii)*0.5 0 0 theta0 dtheta0 0 dphi(ii)]; 
    [Time, X] = ode45(@DerivativesL, [0 15], IC);
    Amp(ii) = (max(abs(X(:,7))) - abs(theta0))/(theta0);
%     if ~isempty(Ie)
%         break
%     end
end
% plot(Time,X(:,7))
% xlabel('Time [sec]'); ylabel('\theta [rad]');
plot(dphi, Amp)
xlabel('$\dot{\phi}$ [rad/sec]','interpreter','Latex'); ylabel('Amplitude','interpreter','latex');
    

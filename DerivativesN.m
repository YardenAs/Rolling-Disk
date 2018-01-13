function [Xdot] = DerivativesL(t, X) %#ok
% The dynamic equations of the rolling disk - Using Newton method

x = X(1); y = X(3); psi = X(5); th = X(7); phi = X(9);
dx = X(2); dy = X(4); dpsi = X(6); dth = X(8); dphi = X(10);
R = 0.5; g = 9.81;
Xdot = [X(2);
        -(dpsi*dy*cos(th) + 2*dpsi*dy*cos(psi)^2*cos(th) + 6*R*dphi*dth*sin(psi)*sin(th) - 5*R*dpsi*dth*cos(th)^2*sin(psi) - 2*dpsi*dx*cos(psi)*cos(th)*sin(psi))/(3*cos(th));
        X(4);
        -(2*dpsi*dx*cos(psi)^2*cos(th) - 3*dpsi*dx*cos(th) - 6*R*dphi*dth*cos(psi)*sin(th) + 5*R*dpsi*dth*cos(psi)*cos(th)^2 + 2*dpsi*dy*cos(psi)*cos(th)*sin(psi))/(3*cos(th));
        X(6);
        -(2*dphi*dth)/cos(th);
        X(8);
        (4*g*sin(th) + 5*R*dpsi^2*cos(th)*sin(th) + 2*R*dphi*dpsi*cos(th) + 4*dpsi*dy*cos(psi)*cos(th) - 4*dpsi*dx*cos(th)*sin(psi))/(5*R);
        X(10);
        -(5*R*dpsi*dth*cos(th)^2 - 6*R*dphi*dth*sin(th) + 2*dpsi*dx*cos(psi)*cos(th) + 2*dpsi*dy*cos(th)*sin(psi))/(3*R*cos(th))];
end





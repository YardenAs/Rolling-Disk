clear all

syms m g Ic R  real                              % problem's parameters
syms x y psi th phi dx dy dpsi dth dphi real     % problem's generalized coordinates
syms ddx ddy ddth ddphi ddpsi ddth real
syms Ft Fn Fs real

Ic = m*R^2/4*[2 0 0; 0 1 0; 0 0 1]; 

q = [x y psi th phi].'; dq = [dx dy dpsi dth dphi].';
ddq = [ddx ddy ddpsi ddth ddphi].';

% Coordinate systems definition
e1p = [cos(psi), sin(psi), 0]; e2p = [-sin(psi), cos(psi), 0]; e3p = [0 0 1];
e1pp = cos(th)*e1p - sin(th)*e3p; e2pp = e2p;
e3pp = sin(th)*e1p + cos(th)*e3p;

% Positions and velocities
rp = [x y 0];
rc = [x y 0] + R*e3pp;
vc = (jacobian(rc,q)*dq).';
vp = (jacobian(rp,q)*dq).';
rcprime = [x*cos(psi) + y*sin(psi) + R*sin(th),...
    y*cos(psi) - x*sin(psi), R*cos(th)];
vcprime = (jacobian(rcprime,q)*dq).' + cross([0 0 dpsi], rcprime);

wpp = [-dpsi*sin(th) - dphi, dth, dpsi*cos(th)].'; % in the body attached coordinate system
wp  = [-dphi*cos(th), dth, dpsi + dphi*sin(th)].';
w = dpsi*e3p + dth*e2p - dphi*e1pp;

% Linear Momentum

acprime = (jacobian(vcprime,q)*dq).' + (jacobian(vcprime,dq)*ddq).' + ...
    cross([0 0 dpsi],vcprime);
acprime = simplify(acprime);

LinMom = [-Fs == m*acprime(1); -Ft == m*acprime(2); Fn - m*g == m*acprime(3)];

% Angular Momentum about COM

Hc = Ic*wpp;
dHc = (jacobian(Hc,q)*dq) + (jacobian(Hc,dq)*ddq) + ...
    cross([-dpsi*sin(th) dth dpsi*cos(th)],Hc).';
dHc = simplify(dHc);

Mc = cross([0 0 -R],[-Fs*cos(th) - Fn*sin(th), -Ft, Fn*cos(th) - Fs*sin(th)]);

AngMom = [dHc(1) == Mc(1); dHc(2) == Mc(2); dHc(3) == Mc(3)];

% No slip constraint

con = simplify(vcprime.' + cross(wp,-[R*sin(th), 0, R*cos(th)]).');
dcon = ((jacobian(con(1:2),q)*dq) + (jacobian(con(1:2),dq)*ddq)) == 0;

% Dynamic Eqautions

sol = solve([AngMom; LinMom; dcon(1:2)], [ddx ddy ddpsi ddth ddphi Ft Fs Fn]);
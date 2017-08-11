syms m g Ic R  real                              % problem's parameters
syms x y psi th phi dx dy dpsi dth dphi real     % problem's generalized coordinates
syms ddx ddy ddth ddphi ddpsi ddth real
syms Ft Fn Fs real
Ic = m*R^2/4*[2 0 0; 0 1 0; 0 0 1]; 

q = [x y psi th phi].'; dq = [dx dy dpsi dth dphi].';
ddq = [ddx ddy ddpsi ddth ddphi].'; LAMBDA = [Ft Fs Fn].';

% Coordinate systems definition
e1p = [cos(psi), sin(psi), 0]; e2p = [-sin(psi), cos(psi), 0]; e3p = [0 0 1];
e1pp = cos(th)*e1p - sin(th)*e3p; e2pp = e2p;
e3pp = sin(th)*e1p + cos(th)*e3p;

% Positions and velocities

rp = [x y 0];
rc = [x y 0] + R*e3pp;
vc = (jacobian(rc,q)*dq).';
vp = (jacobian(rp,q)*dq).';
vcp = [dx*cos(psi) + dy*sin(psi) + R*dth*cos(th), dy*cos(psi) - dx*sin(psi)...
    + dpsi*R*sin(th), -R*dth*sin(th)]; % velocity of COM in ei' coordinate system

wpp = [-dpsi*sin(th) + dphi*cos(th), dth, dpsi*cos(th) - dphi*sin(th)]; % in the body attached coordinate system
wp  = [-dphi*cos(th), dth, dpsi];

% Kinetic energy

T = 1/2*m*(vc*vc.') + 1/2*wpp*Ic*wpp.';

% Potential energy 

P = -[0 0 -m*g]*rc.';

% Non - holonomic constraints (on ei' coordinate system)

con = simplify((vcp + cross(wp,[R*sin(th), 0, R*cos(th)])).');

% Equation's matrices

W = simplify(jacobian(con,dq));
M = simplify(hessian(T,dq));
G = simplify(jacobian(P,q).');

C  = sym(zeros(length(q),1));

for ii = 1:length(q)
    sub(ii) = diff(T,dq(ii)); %#ok
    for jj = 1:4
        subb = diff(sub(ii),q(jj));
        C(ii) = C(ii) + subb*dq(jj); 
    end
    C(ii) = C(ii) - diff(T,q(ii));
end
C = simplify(C);                  % C vector



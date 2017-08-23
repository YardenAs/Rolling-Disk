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
vcprime = [dx*cos(psi) + dy*sin(psi) + R*dth*cos(th), dy*cos(psi) - dx*sin(psi)...
    + dpsi*R*sin(th), -R*dth*sin(th)]; % COM velocity in ei' coordinate system

wpp = [-dpsi*sin(th) - dphi, dth, dpsi*cos(th)].'; % in the body attached coordinate system
wp  = [-dphi*cos(th), dth, dpsi + dphi*sin(th)].';
w = dpsi*e3p + dth*e2p - dphi*e1pp;

% Kinetic energy
T = 1/2*m*(vc*vc.') + 1/2*wpp.'*Ic*wpp;

% Potential energy 
V = -[0 0 -m*g]*rc.';

% Non - holonomic constraints (on ei' coordinate system)
conT = simplify((vcprime + cross(wp,-[R*sin(th), 0, R*cos(th)])).');
constraint = conT(1:2);

% Equation's matrices
W = simplify(jacobian(constraint,dq));
M = simplify(hessian(T,dq));
G = simplify(jacobian(V,q).');

C  = sym(zeros(length(q),1));
for ii = 1:length(q)
    sub(ii) = diff(T,dq(ii)); %#ok
    for jj = 1:length(q)
        subb = diff(sub(ii),q(jj));
        C(ii) = C(ii) + subb*dq(jj); 
    end
    C(ii) = C(ii) - diff(T,q(ii));
end
C = simplify(C); % C vector

n = size(W);
Wdot = sym(zeros(n));
for i = 1:n(1)
    Wdot(i,:) = (jacobian(W(i,:),q)*dq).';
end
Wdot = simplify(Wdot);

a = [M -W.';
    W zeros(2,2)];
    
b = [-C - G;
    -Wdot*dq];

sol = simplify(a\b);

% Substitute the non holonomic constraints into ddth and ddpsi
solCon = solve(conT == 0, [dx dy]);
ddth = subs(sol(4),[dx dy], [solCon.dx solCon.dy]);
ddpsi = subs(sol(3),[dx dy], [solCon.dx solCon.dy]);

%% Linarization %% 
syms dphi0
ddpsi = subs(ddpsi,[dphi dth th],[dphi0 0 0]) +...
    subs(jacobian(ddpsi,[dphi, dth th]),[dphi dth th], [dphi0 0 0])*[dphi dth th].';
ddth = subs(ddth,[dphi dpsi th],[dphi0 0 0]) +...
    subs(jacobian(ddth,[dphi, dpsi th]),[dphi dpsi th], [dphi0 0 0])*[dphi dpsi th].';
dddth = jacobian(ddth,q)*dq + jacobian(ddth,dq)*ddq;
syms m g Ic R  real                              % problem's parameters
syms x y psi th phi dx dy dpsi dth dphi real     % problem's generalized coordinates
syms ddx ddy ddth ddphi ddpsi ddth real
syms alpha             

Ic = alpha*m*R^2/2*[2 0 0; 0 1 0; 0 0 1]; 

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
solCon = solve(conT == 0, [dx dy]);

% Substitute the non holonomic constraint
T = subs(T,[dx dy], [solCon.dx solCon.dy]);
V = subs(V,[dx dy], [solCon.dx solCon.dy]);
q(1:2) = [];
dq(1:2) = [];
ddq(1:2) = [];

% Equation's matrices
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

sol = M\(- C - G);
%% Derive ddth = ddth(th, C1, C2) formula %%
syms C1 C2 real
L = T - V;
ddth = sol(2);
ddpsi = sol(3);
% ConservationMom1 = diff(L,dpsi) == C1;
% ConservationMom2 = diff(L,dphi) == C2;
% solConMom = solve([ConservationMom1 ConservationMom2], [dpsi, dphi]);
% ddth = subs(ddth, [dpsi dphi], [solConMom.dpsi solConMom.dphi]);
% 
% %% ddth when the disk is upright %%
% syms dphi0
% C10 = subs(diff(L,dpsi),[th dphi dpsi],[0 dphi0 0]);
% C20 = subs(diff(L,dphi),[th dphi dpsi],[0 dphi0 0]);
% ddth = subs(ddth,[C1 C2], [C10 C20]);
% % linearize about th = 0
% ddth = subs(ddth,th,0) + th*subs(diff(ddth,th),th,0);

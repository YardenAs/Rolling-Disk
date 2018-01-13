function [value, isterminal, direction] = isStable(~, X)
% rolling disk - instability event if the amplitude of theta is twice larger
% than the initial value - the disk is unstable
global theta0
global dtheta0
isterminal = ones(1,2);
direction  = -ones(1,2);
if theta0 ~= 0
%     value(1) = 1 - (abs(X(7)) - abs(theta0))/abs(theta0);
value(1) = 2*abs(theta0) - abs(X(7));
else value(1) = 1;
end
if dtheta0 ~= 0
    value(2) = 2*abs(dtheta0) - abs(X(8));
else value(2) = 1;
end
end


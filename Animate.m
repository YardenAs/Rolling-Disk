function Animate(Time ,X)
% This function animates a rolling disk

x = X(:,1); y = X(:,3); psi = X(:,5); th = X(:,7); phi = X(:,9);
dpsi = X(:,6); dth = X(:,8); dphi = X(:,10); %#ok
R = 0.5;
ang = linspace(0, 2*pi, 40).';
circle = [R*cos(ang) R*sin(ang)];
len = length(X);

% Coordinate systems definition
e1p = [cos(psi), sin(psi), zeros(len,1)];
e2p = [-sin(psi), cos(psi), zeros(len,1)];
e3p = [zeros(len,1) zeros(len,1) ones(len,1)];  
e2pp = e2p;
e3pp = zeros(len,3); e1pp = zeros(len,3); c = zeros(len,3);
xcircle = zeros(len,40); ycircle = zeros(len,40); zcircle = zeros(len,40);
pcircle = [zeros(40,1) circle(:,1) circle(:,2)];
for ii = 1:len
    e3pp(ii,:) = e1p(ii,:)*sin(th(ii)) + cos(th(ii))*e3p(ii,:); 
    e1pp(ii,:) = e1p(ii,:)*cos(th(ii)) - sin(th(ii))*e3p(ii,:);
    c(ii,:) = [x(ii) y(ii) 0] + R*e3pp(ii,:);
    xcircle(ii,:) = c(ii,1) + circle(:,1)*e2pp(ii,1) + circle(:,2)*e3pp(ii,1); 
    ycircle(ii,:) = c(ii,2) + circle(:,1)*e2pp(ii,2) + circle(:,2)*e3pp(ii,2);
    zcircle(ii,:) = c(ii,3) + circle(:,1)*e2pp(ii,3) + circle(:,2)*e3pp(ii,3);
    e3ppp(ii,:) = sin(phi(ii))*e2pp(ii,:) + cos(phi(ii))*e3pp(ii,:); %#ok
    a(ii,:) = c(ii,:) + R*e3ppp(ii,:); %#ok
end

figure
set(gcf, 'color', 'w')
view(3);
grid on
path = animatedline('color', 'b', 'linewidth', 2);
disk = patch('xdata', pcircle(:,1), 'ydata', pcircle(:,2), 'zdata', pcircle(:,3), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);
pointA = line('xdata', a(1,1), 'ydata', a(1,2), 'zdata', a(1,3)-R, 'marker', 'o', 'color', 'r', 'markerfacecolor', 'r', 'linewidth', 1);

axis equal
xlim([min(x) - R, max(x) + R]); ylim([min(y) - R, max(y) + R]); zlim([-1.2*2*R 1.2*2*R]);

for ii = 1:len
    
    set(disk,'xdata',xcircle(ii,:), 'ydata', ycircle(ii,:), 'zdata', zcircle(ii,:));
    addpoints(path, x(ii), y(ii), 0);
    set(pointA, 'xdata', a(ii,1), 'ydata', a(ii,2), 'zdata', a(ii,3));
    drawnow
%     pause(0.1);
end



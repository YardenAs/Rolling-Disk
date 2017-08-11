function Animate(X)
% This function animates a rolling disk

x = X(:,1); y = X(:,2); psi = X(:,3); th = X(:,4); phi = X(:,5);
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
for ii = 1:len
    e3pp(ii,:) = e1p(ii,:)*sin(th(ii)) + cos(th(ii))*e3p(ii,:); 
    e1pp(ii,:) = e1p(ii,:)*cos(th(ii)) - sin(th(ii))*e3p(ii,:);
    c(ii) = [x y zeros(len,1)] + R*e3pp;
    xcircle(ii,:) = c(ii,1) + circle(:,1)*e2pp(ii,1) + circle(:,2)*e3pp(ii,1); 
    ycircle(ii,:) = c(ii,2) + circle(:,1)*e2pp(ii,2) + circle(:,2)*e3pp(ii,2);
    zcircle(ii,:) = c(ii,3) + circle(:,1)*e2pp(ii,3) + circle(:,2)*e3pp(ii,3);
end

figure
set(gcf, 'color', 'w')
plot3(x1(1), x2(1), x3(1));
xlabel('\itx\rm_{1} (m)')
set(gca, 'xdir', 'reverse')
ylabel('\itx\rm_{2} (m)')
set(gca, 'ydir', 'reverse')
zlabel('\itx\rm_{3} (m)            ', 'rotation', 0)
axis equal
xlim([min(x1)-r, max(x1)+r])
ylim([min(x2)-r, max(x2)+r])
zlim([0, 1.2*(2*r)])
grid on

disk = patch('xdata', xcircle(1,:), 'ydata', ycircle(1,:), 'zdata', zcircle(1,:), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);


for ii = 1:len
    set(disk, 'xdata', xcircle(:,ii), 'ydata', ycircle(:,ii), 'zdata', zcircle(:,ii));
    drawnow   
end 



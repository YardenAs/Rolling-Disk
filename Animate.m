function Animate(X)
% This function animates a rolling disk

x = X(:,1); y = X(:,3); psi = X(:,5); th = X(:,7); phi = X(:,9);
dpsi = X(:,6); dth = X(:,8); dphi = X(:,10); %#ok
R = 0.5;
ang = linspace(0, 2*pi, 40).';
circle = [R*cos(ang) R*sin(ang)];
pcircle = [zeros(40,1) circle(:,1) circle(:,2)];

% Coordinate systems definition
len = length(X);
e1p = [cos(psi), sin(psi), zeros(len,1)];
e2p = [-sin(psi), cos(psi), zeros(len,1)];
e3p = [zeros(len,1) zeros(len,1) ones(len,1)];  
e3pp = zeros(len,3); e1pp = zeros(len,3); c = zeros(len,3); e2pp = e2p; %#ok

% Calculate ei_pp and COM for each timestep.
for ii = 1:len
    e3pp(ii,:) = e1p(ii,:)*sin(th(ii)) + cos(th(ii))*e3p(ii,:); 
    e1pp(ii,:) = e1p(ii,:)*cos(th(ii)) - sin(th(ii))*e3p(ii,:);
    c(ii,:) = [x(ii) y(ii) 0] + R*e3pp(ii,:);
end


figure;
set(gcf, 'color', 'w');
view(3);
grid on
hg = hgtransform('Parent',gca);
path = animatedline('color', 'b', 'linewidth', 2);
disk(1) = patch('xdata', pcircle(:,1), 'ydata', pcircle(:,2),...
    'zdata', pcircle(:,3), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);
disk(2) = line('xdata', pcircle(1,1), 'ydata', pcircle(1,2), 'zdata',...
    pcircle(1,3), 'marker', 'o', 'color', 'r', 'markerfacecolor', 'r', 'linewidth', 1);
set(disk, 'Parent', hg);
axis equal

xlim([min(x) - 2*R, max(x) + 2*R]); ylim([min(y) - 2*R, max(y) + 2*R]);
zlim([0 1.2*2*R]);

for ii = 1:len
    T = makehgtform('translate', c(ii,:));
    R = makehgtform('axisrotate', e1pp(ii,:), -phi(ii), 'axisrotate',...
        e2p(ii,:), th(ii), 'zrotate', psi(ii));
    % Important note about R: note that the first argument is e1pp then e2p
    % then e3! e1pp is defined by ei_p and ei_p is defined by ei. if we
    % change the order of the arguments, we will get rubish.
    set(hg, 'Matrix', T*R); % T times R and not R times T!
    addpoints(path, x(ii), y(ii), 0);
    drawnow
%     pause(0.01)
end




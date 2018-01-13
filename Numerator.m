function [Mone] = Numerator(dphi0, ths)
%Numerator evaluates the numerator of ddeta
g = 9.81; R = 0.5; 
Mone = 10*R*g - 4*dphi0.*(R^3*cos(ths)*(20*g*cos(ths).^2 - 20*g + 9*R*dphi0.^2*cos(ths))).^(1/2) + 9*R^2*dphi0.^2.*cos(ths).^3 - 10*R*g*cos(ths).^4 - 18*R^2*dphi0.^2.*cos(ths) + 7*dphi0.*cos(ths).^2.*(R^3*cos(ths).*(20*g*cos(ths).^2 - 20*g + 9*R*dphi0.^2.*cos(ths))).^(1/2);
end


function [r] = Multilateration3D(v1,r1,v2,r2,v3,r3)
%3-D MULTILATERATION Performs 3-dimensional multilateration
%   Takes as inputs [x; y; z] location and range-to-target for three  
%   detectors, gives as output estimated [x; y; z] of target.

% Translate v1 to (0,0,0)
w2 = v2-v1;
w3 = v3-v1;

% Rotate v2 onto x-axis
a = sqrt(sum(w2.^2));
b = sqrt(w2(1)^2 + w2(3)^2);

if b == 0
    Ry = [0 0 1; 0 1 0; 1 0 0];
else
    Ry = [w2(1)/b, 0, w2(3)/b; 0, 1, 0; -w2(3)/b, 0, w2(1)/b];
end

Rz = [b/a, w2(2)/a, 0; -w2(2)/a, b/a, 0; 0, 0, 1];

R1 = Rz*Ry;

u2 = R1*w2;
u3 = R1*w3;

% Rotate v3 into x-y plane
c = sqrt(u3(2)^2 + u3(3)^2);
if c == 0
    R2 = eye(3);
else
    R2 = [1, 0, 0; 0, u3(2)/c, u3(3)/c; 0, -u3(3)/c, u3(2)/c];
end

t2 = R2*u2;
t3 = R2*u3;

% Calculate target location for new coordinates
x = (r1^2 - r2^2 + t2(1)^2)/(2*t2(1));
y = (r1^2 - r3^2 + t3(1)^2 - 2*t3(1)*x + t3(2)^2)/(2*t3(2));
z = sqrt(abs(r1^2 - x^2 - y^2));

r = [x;y;z];

% Inverse transformation to original coordinate system
r = R1\(R2\r) + v1;

r(3) = abs(r(3));

r(isnan(r)) = 0;

end

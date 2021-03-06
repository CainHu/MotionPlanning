function [pos, vel, acc, jerk, snap] = bezier5(u, P)

n = 5;
P0 = P(1);
P1 = P(2);
P2 = P(3);
P3 = P(4);
P4 = P(5);
P5 = P(6);

% a = 1-u;
% a2 = a * a;
% a3 = a2 * a;
% a4 = a3 * a;
% a5 = a4 * a;
% b = u;
% b2 = b * b;
% b3 = b2 * b;
% b4 = b3 * b;
% b5 = b4 * b;

pos = [(1-u)^5, 5*u*(1-u)^4, 10*u^2*(1-u)^3, 10*u^3*(1-u)^2, 5*u^4*(1-u), u^5] * [P0, P1, P2, P3, P4, P5]';
vel = [(1-u)^4, 4*u*(1-u)^3, 6*u^2*(1-u)^2, 4*u^3*(1-u), 1*u^4] * [P1-P0, P2-P1, P3-P2, P4-P3, P5-P4]' * n;
acc = [(1-u)^3, 3*u*(1-u)^2, 3*u^2*(1-u), u^3] * [P2-2*P1+P0, P3-2*P2+P1, P4-2*P3+P2, P5-2*P4+P3]' * n*(n-1);
jerk = [(1-u)^2, 2*u*(1-u), u^2] * [P3-3*P2+3*P1-P0, P4-3*P3+3*P2-P1, P5-3*P4+3*P3-P2]' * n*(n-2)*(n-2);
snap = [(1-u), u] * [P4-4*P3+6*P2-4*P1+P0, P5-4*P4+6*P3-4*P2+P1]' * n*(n-2)*(n-2)*(n-3);


end


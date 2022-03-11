syms c1 c2 c3 c4 c5 y1 y2 y3 sigma_ini p_n p_n_1 f_o a b dt k L
eqn1 = c1 * (-y1-sign(y3+c3) * (sigma_ini-p_n)/(k*L)*(f_o + a * log(y3+c3) + b * log(y2)));
eqn2 = c4 - abs(y3+c3) * y2 + c5 * y2 * (p_n-p_n_1)/dt;
eqn3 = y3;

deqn1 = [diff(eqn1,y1) diff(eqn1,y2) diff(eqn1,y3)];
deqn2 = [diff(eqn2,y1) diff(eqn2,y2) diff(eqn2,y3)];
deqn3 = [diff(eqn3,y1) diff(eqn3,y2) diff(eqn3,y3)];

pretty(deqn1)
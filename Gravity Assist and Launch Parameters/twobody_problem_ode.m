function ds = twobody_problem_ode( t, s, muP)

r= s(1:3);
v= s(4:6);

r_norm=norm(r);

ds=[v(1);v(2);v(3);-muP/r_norm^3*r(1);-muP/r_norm^3*r(2);-muP/r_norm^3*r(3)];

end
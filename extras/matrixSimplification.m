syms argPeriapsis trueAnomaly inclination RAAN

A1 = [cos(argPeriapsis + trueAnomaly), sin(argPeriapsis + trueAnomaly), 0; -sin(argPeriapsis + trueAnomaly), cos(argPeriapsis + trueAnomaly), 0; 0, 0, 1]
A2 = [1, 0, 0; 0, cos(inclination), sin(inclination); 0, -sin(inclination), cos(inclination)]
A3 = [cos(RAAN), sin(RAAN), 0; -sin(RAAN), cos(RAAN), 0; 0, 0, 1]

inv(mtimes(A1, mtimes(A2, A3)))

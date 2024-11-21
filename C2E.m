function [orb] = C2E(rv, mu)
    % Unpack the Cartesian state vector
    r_vec = rv(1:3); % Position vector (km)
    v_vec = rv(4:6); % Velocity vector (km/s)
    
    % Magnitudes of position and velocity
    r = norm(r_vec);
    v = norm(v_vec);
    
    % Specific angular momentum
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    
    % Inclination
    i = acos(h_vec(3) / h);
    
    % Node line
    N_vec = cross([0; 0; 1], h_vec);
    N = norm(N_vec);
    
    % Right Ascension of Ascending Node (RAAN)
    if N_vec(2) >= 0
        Omega = acos(N_vec(1) / N);
    else
        Omega = 2 * pi - acos(N_vec(1) / N);
    end
    
    % Eccentricity vector
    e_vec = (1 / mu) * ((v^2 - mu / r) * r_vec - dot(r_vec, v_vec) * v_vec);
    e = norm(e_vec);
    
    % Argument of periapsis
    if e_vec(3) >= 0
        omega = acos(dot(N_vec, e_vec) / (N * e));
    else
        omega = 2 * pi - acos(dot(N_vec, e_vec) / (N * e));
    end
    
    % True anomaly
    if dot(r_vec, v_vec) >= 0
        nu = acos(dot(e_vec, r_vec) / (e * r));
    else
        nu = 2 * pi - acos(dot(e_vec, r_vec) / (e * r));
    end
    
    % Semi-major axis
    a = 1 / (2 / r - v^2 / mu);
    
    % Mean anomaly (M) approximation (from eccentric anomaly)
    E = atan2(sqrt(1 - e^2) * sin(nu), e + cos(nu));
    M = E - e * sin(E);
    
    % Output the orbital elements
    orb = [a; e; i; Omega; omega; M];
end

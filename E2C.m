function [rv] = E2C(orb, mu)
    % Unpack the orbital elements
    a = orb(1); % Semi-major axis (km)
    e = orb(2); % Eccentricity
    i = orb(3); % Inclination (rad)
    Omega = orb(4); % Right Ascension of Ascending Node (rad)
    omega = orb(5); % Argument of periapsis (rad)
    M = orb(6); % Mean anomaly (rad)
    
    % Solve Kepler's equation for eccentric anomaly
    E = M; % Initial guess for eccentric anomaly
    for iter = 1:100
        E_new = M + e * sin(E);
        if abs(E_new - E) < 1e-10
            break;
        end
        E = E_new;
    end
    
    % Compute the true anomaly
    nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    
    % Distance
    r = a * (1 - e * cos(E));
    
    % Position in the perifocal frame
    r_perifocal = [r * cos(nu); r * sin(nu); 0];
    
    % Velocity in the perifocal frame
    h = sqrt(mu * a * (1 - e^2)); % Specific angular momentum
    v_perifocal = [-mu/h * sin(nu); mu/h * (e + cos(nu)); 0];
    
    % Rotation matrix from perifocal to geocentric equatorial frame
    R3_Omega = [cos(Omega), sin(Omega), 0;
               -sin(Omega), cos(Omega), 0;
                0, 0, 1];
    R1_i = [1, 0, 0;
            0, cos(i), sin(i);
            0, -sin(i), cos(i)];
    R3_omega = [cos(omega), sin(omega), 0;
               -sin(omega), cos(omega), 0;
                0, 0, 1];
    
    Q_pX = (R3_Omega') * (R1_i') * (R3_omega'); % Final rotation matrix
    
    % Position and velocity in geocentric equatorial frame
    r_xyz = Q_pX * r_perifocal;
    v_xyz = Q_pX * v_perifocal;
    
    rv = [r_xyz; v_xyz]; % Output position and velocity
end

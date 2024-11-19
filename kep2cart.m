function [r, v] = kep2cart(kep, mu)
    % Convert Keplerian elements to Cartesian coordinates
    % INPUTS:
    %   kep - [a e i Om om Mo], Keplerian elements
    %   mu  - Gravitational parameter of the central body
    % OUTPUTS:
    %   r   - Position vector [km]
    %   v   - Velocity vector [km/s]

    % Extract elements
    a = kep(1); % Semi-major axis
    e = kep(2); % Eccentricity
    i = kep(3); % Inclination
    Om = kep(4); % Longitude of ascending node
    om = kep(5); % Argument of periapsis
    Mo = kep(6); % Mean anomaly at epoch

    % Solve Kepler's Equation for Eccentric Anomaly using Newton-Raphson method
    E = Mo; % Initial guess for Eccentric Anomaly
    tol = 1e-8; % Convergence tolerance
    max_iter = 100; % Maximum number of iterations

    for iter = 1:max_iter
        f = E - e * sin(E) - Mo; % Kepler's equation
        f_prime = 1 - e * cos(E); % Derivative of Kepler's equation
        E_new = E - f / f_prime; % Update using Newton-Raphson method
        if abs(E_new - E) < tol
            E = E_new;
            break; % Converged
        end
        E = E_new;
    end

    % Compute true anomaly (nu)
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Compute distance (r_mag)
    r_mag = a * (1 - e * cos(E));

    % Position in the orbital plane (perifocal coordinates)
    r_perifocal = [r_mag * cos(nu); r_mag * sin(nu); 0];

    % Velocity in the orbital plane (perifocal coordinates)
    h = sqrt(mu * a * (1 - e^2)); % Specific angular momentum
    v_perifocal = (mu / h) * [-sin(nu); e + cos(nu); 0];

    % Rotation matrices to convert from perifocal to inertial frame
    R3_W = [cos(Om), -sin(Om), 0; sin(Om), cos(Om), 0; 0, 0, 1]; % Rotation about z by Om
    R1_i = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)]; % Rotation about x by i
    R3_w = [cos(om), -sin(om), 0; sin(om), cos(om), 0; 0, 0, 1]; % Rotation about z by om

    % Combined rotation matrix
    Q_pX = R3_W * R1_i * R3_w;

    % Position and velocity in the inertial frame
    r = Q_pX * r_perifocal; % Position in inertial frame
    v = Q_pX * v_perifocal; % Velocity in inertial frame
end

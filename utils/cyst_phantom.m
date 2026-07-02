function [positions, amp] = cyst_phantom (N)

    x_size = 50/1000;   %  Width of phantom [mm]
    y_size = 10/1000;   %  Transverse width of phantom [mm]
    z_size = 60/1000;   %  Height of phantom [mm]
    z_start = 30/1000;  %  Start of phantom surface [mm];

    % Create the general scatterers
    x = (rand(N, 1) - 0.5) * x_size;
    y = (rand(N, 1) - 0.5) * y_size;
    z = rand(N, 1) * z_size + z_start;

    % Generate the amplitudes with a Gaussian distribution
    amp = randn(N,1);

    % Make the cyst and set the amplitudes to zero inside
    r = 10/1000;     % Radius of cyst [m]
    xc = 0/1000;    % Place of cyst [m]
    zc = 25/1000 + z_start;

    inside = (((x-xc).^2 + (z-zc).^2) < r^2);
    amp = amp .* (1-inside);

    % Place the point scatterers in the phantom
    for i=N-5:N
        x(i) = -15/1000;
        y(i) = 0;
        z(i) = z_start + (10+5*10)/1000 + (i-N)*10/1000;
        amp(i) = 20;
    end
    % Return the variables
    positions=[x y z];
end

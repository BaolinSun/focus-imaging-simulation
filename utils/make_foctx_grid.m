%% Generate a focused pixel grid based on input parameters
function grid = make_foctx_grid(rlims, dr, dirs)

    r = rlims(1) : dr : rlims(2);
    t = dirs(:, 1);
    [tt, rr] = meshgrid(t, r);
    rr = rr';
    tt = tt';

    xx = rr .* sin(tt);
    zz = rr .* cos(tt);
    yy = zeros(size(xx));

    grid = cat(3, xx, yy, zz);
end
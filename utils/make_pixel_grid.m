%% Generate a Cartesian pixel grid based on input parameters.
function grid = make_pixel_grid(xlims, zlims, dx, dz)

    x = xlims(1):dx:xlims(2);
    z = zlims(1):dz:zlims(2);

    [xx, zz] = meshgrid(x, z);
    yy = zeros(size(xx));

    grid = cat(3, xx, yy, zz);
end
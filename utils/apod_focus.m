%% Compute rect apodization to user-defined pixels for desired f-number
function apod = apod_focus(grid, ele_pos, fnum, hamming_win)

    min_width = 1e-3;

    ppos = reshape(grid, [1, size(grid, 2), 3]);
    epos = reshape(ele_pos, [size(ele_pos, 1), 1, 3]);

    v = ppos - epos;
    v_x = v(:, :, 1);
    v_z = v(:, :, 3);

    mask_part1 = abs(v_z ./ v_x) >= fnum;    % 动态孔径条件
    mask_part2 = abs(v_x) <= min_width;      % 最小孔径条件

    mask = mask_part1 | mask_part2;

    win = repmat(hamming_win, 1, size(grid, 2));

    apod = mask .* win;
end
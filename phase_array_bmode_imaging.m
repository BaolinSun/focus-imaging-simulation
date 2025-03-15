% phase_array_bmode_imaging differs from phase_array_focus_imaging only in the use of actual rf data

% Field II Phase Array B-mode
clear all; close all;

%% Init Field II
path(path, 'D:\MyProjects\matlab\Field_II_ver_3_30_windows');
field_init(-1);

%% 子函数路径
addpath('utils')
addpath("probe")

%% 相控阵探头参数
probe = Probe('phase array');
f0 = 2.5e6;              % 中心频率 3.5 MHz
element_num = 64;        % 阵元数量
height = 5e-3;           % 阵元高度 (m)
width = 0.12e-3;         % 阵元宽度 (m)
kerf = 0.18e-3;          % 阵元间距 (m)
pitch = width + kerf;
focus = 60e-3;           % 发射聚焦深度 (m)
c = 1540;                % 声速 (m/s)
fs = 20e6;              % 采样频率 (Hz)
Ts = 1/fs;               % 采样间隔 (s)

% 阵元参数
x_ele = ([0:element_num-1]-(element_num-1)/2).*pitch;
z_ele  = zeros(1,length(x_ele));
probe.N_elements = element_num;
probe.pitch = pitch;
probe.x_ele = x_ele;
probe.y_ele = z_ele;
probe.z_ele = z_ele;
probe.ele_pos = [x_ele; z_ele; z_ele]';


%% 扇扫参数设置
F = 120e-3;
% scan_angle_deg = -30:2:30;
angles = linspace(-45, 45, 64);
theta = deg2rad(angles);
num_line = length(angles);
raw_data = cell(1, num_line); % 存储每条扫描线的多阵元原始数据
tstart = zeros(num_line);


%% 读取csv文件
for i = 1:num_line
    rfdata = readmatrix(['rawdata\rfdata\rfdata_1_', num2str(i), '.csv']);
    % rfdata = (rfdata - 512) / 512;
    rfdata = bandpass_filter(rfdata);
    raw_data{i} = rfdata;
    tstart(i) = 0;
end


%% 接收波束合成
rx_num_line = 64;
parallel_beam = rx_num_line / num_line;

ele_pos = probe.ele_pos;

oris = zeros(rx_num_line, 3);
oris(:, 1) = linspace(0, rx_num_line - 1, rx_num_line) * pitch * 64.0 / rx_num_line;
oris(:, 1) = oris(:, 1) - mean(oris(:, 1));
tx_ori = oris;

rx_dir = [theta; zeros(size(theta))]';

rmax = 158e-3;
rlims = [0, rmax];
wvln = c / f0;
dr = wvln / 4;
r = rlims(1) : dr : rlims(2);
grid = make_foctx_grid(rlims, dr, rx_dir);
grid_s = size(grid);
nx = grid_s(1);
nz = grid_s(2);
das = zeros(nx, nz);
foc = zeros(rx_num_line, nz);

hann_window = hann(element_num);

for i = 1:rx_num_line
    data_line = ceil(i / parallel_beam);
    data = raw_data{i}';
    % data = data .* hann_window;

    txdel = vecnorm(squeeze(grid(i, :, :)) - squeeze(tx_ori(i, :, :)), 2, 2)';
    rxdel = sqrt(sum((reshape(grid(i, :, :), [], 1, 3) - reshape(ele_pos, [1, size(ele_pos)])).^2, 3))';
    delays = ((txdel + rxdel) / c - tstart(data_line)) * fs;

    xc = 1 : size(data, 2);
    for j = 1 : element_num
        analytic_signal = hilbert(data(j, :));
        foc(j, :) = interp1(xc, analytic_signal, delays(j, :), 'linear', 0.0);
    end

    apods = apod_focus(grid(i, :, :), ele_pos, 1);
    foc = foc .* apods;

    das(i, :) = sum(foc);
end

%% Scan convert
xlims = rlims(2) * [-0.7, 0.7];
zlims = rlims(2) * [0, 1];
img_grid = make_pixel_grid(xlims, zlims, wvln / 2, wvln / 2);

grid_x = grid(:, :, 1);
grid_x = grid_x(:);
grid_y = grid(:, :, 3);
grid_y = grid_y(:);
img_grid_x = img_grid(:, :, 1);
img_grid_x = img_grid_x(:);
img_grid_y = img_grid(:, :, 3);
img_grid_y = img_grid_y(:);
bimgsc = griddata(grid_x, grid_y, das(:), img_grid_x, img_grid_y, 'cubic');
% bimgsc(isnan(bimgsc)) = 1e-10;
bimg = reshape(bimgsc, size(img_grid, 1), size(img_grid, 2));


drange = 50;
bimg = abs(bimg);
bimg = 20 * log10(bimg);
bimg = bimg - max(bimg(:));
extent = [img_grid(1, 1, 1), img_grid(1, end, 1), img_grid(1, 1, 3), img_grid(end, 1, 3)] * 1e3;
figure;
imagesc(extent([1 2]), extent([3 4]), bimg);
colormap gray;
caxis([-drange 0]);
% set(gca, 'YDir', 'normal');
xlabel('Lateral distance [mm]');
ylabel('Axis distance [mm]');
axis image;
% colorbar;




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

%% Generate a Cartesian pixel grid based on input parameters.
function grid = make_pixel_grid(xlims, zlims, dx, dz)

    x = xlims(1):dx:xlims(2);
    z = zlims(1):dz:zlims(2);

    [xx, zz] = meshgrid(x, z);
    yy = zeros(size(xx));

    grid = cat(3, xx, yy, zz);
end

%% Compute rect apodization to user-defined pixels for desired f-number
function apod = apod_focus(grid, ele_pos, fnum)

    min_width = 1e-3;

    ppos = reshape(grid, [1, size(grid, 2), 3]);
    epos = reshape(ele_pos, [size(ele_pos, 1), 1, 3]);

    v = ppos - epos;
    v_x = v(:, :, 1);
    v_z = v(:, :, 3);

    mask_part1 = abs(v_z ./ v_x) >= fnum;    % 动态孔径条件
    mask_part2 = abs(v_x) <= min_width;      % 最小孔径条件

    mask = mask_part1 | mask_part2;

    element_num = size(ele_pos, 1);
    hamming_win = hamming(element_num);
    win = repmat(hamming_win, 1, size(grid, 1));

    apod = mask .* win;
end
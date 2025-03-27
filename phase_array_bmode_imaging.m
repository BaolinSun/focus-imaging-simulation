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
tstart = zeros(1, num_line);


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

tx_ori = zeros(rx_num_line, 3);

rx_dir = [theta; zeros(size(theta))]';

rmax = 158e-3;
rlims = [0, rmax];
wvln = c / f0;
dr = wvln / 10;
r = rlims(1) : dr : rlims(2);
grid = make_foctx_grid(rlims, dr, rx_dir);
grid_s = size(grid);
nx = grid_s(1);
nz = grid_s(2);
das = zeros(nx, nz);
foc = zeros(rx_num_line, nz);

hamming_win = hamming(element_num);

segment_length = 256;
cutoff_freq = 1e6;
n = 70;

for i = 1:rx_num_line
    data_line = ceil(i / parallel_beam);
    data = raw_data{i}';
    % data = data .* hann_window;

    % txdel = vecnorm(squeeze(grid(i, :, :)) - squeeze(tx_ori(i, :, :)), 2, 2)';   % (1026x3) - (1x3)
    txdel = sqrt(sum((squeeze(grid(i, :, :)) - squeeze(tx_ori(i, :, :))).^2, 2))';   % (1026x3) - (1x3)
    rxdel = sqrt(sum((reshape(grid(i, :, :), [], 1, 3) - reshape(ele_pos, [1, size(ele_pos)])).^2, 3))';   % (1026x1x3) - (1x64x3)
    delays = ((txdel + rxdel) / c - tstart(data_line)) * fs;

    xc = 1 : size(data, 2);
    for j = 1 : element_num
        % analytic_signal = hilbert(data(j, :));
        foc(j, :) = interp1(xc, data(j, :), delays(j, :), 'linear', 0.0);
    end

    apods = apod_focus(grid(i, :, :), ele_pos, 1, hamming_win);
    foc = foc .* apods;

    beamdata = sum(foc);

    num_segments = ceil(length(beamdata) / segment_length);
    If = zeros(size(beamdata));
    Qf = zeros(size(beamdata));

    for seg = 1:num_segments
        start_idx = (seg-1)*segment_length + 1;
        end_idx = min(seg*segment_length, length(beamdata));

        tdata = beamdata(start_idx:end_idx);
        mag = abs(fft(tdata));
        m = size(mag,2);
        [p k] = max(mag(1:fix(m/2)));
        w = 2 * pi * k/m * fs;
        x = (0:m-1) / fs;
        I = cos(w*x).*tdata;
        Q = sin(w*x).*tdata;        

        % 计算当前段的平均深度
        wn = cutoff_freq / (fs/2);
        fir = fir1(n, wn, 'low', hamming(n+1));

        delay = floor(n/2);
        If_seg = fftfilt(fir, [I, zeros(1, delay)]); % 补零避免截断
        Qf_seg = fftfilt(fir, [Q, zeros(1, delay)]);

        % 存储滤波结果（跳过前delay个点）
        stored_start = start_idx + delay;
        stored_end = min(end_idx + delay, length(If));
        valid_length = stored_end - stored_start + 1;
        If(stored_start:stored_end) = If_seg(delay+1:delay+valid_length);
        Qf(stored_start:stored_end) = Qf_seg(delay+1:delay+valid_length);
    end

    Fdata = sqrt(If.^2 + Qf.^2);
    das(i, :) = Fdata;
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




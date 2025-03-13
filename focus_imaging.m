% Init the Field II simulation system
path(path, 'D:\MyProjects\matlab\Field_II_ver_3_30_windows');
field_init(-1);

wavetype = 'fcous_wave'; %plane_wave, diverging_wave, fcous_wave

%% 探头参数
probe = Probe('linear array');
c = 1540;
f0 = 2500000;
N_elements = 64;
width = 0.12e-3;
height = 5e-3;
kerf = 0.18e-3;
pitch = width+kerf;
focus = 20e-3;
fs = 25e6;
set_sampling(fs)
Th = xdc_linear_array (N_elements, width, height, kerf, 1, 10, [0, 0, focus]);
Rh = xdc_linear_array (N_elements, width, height, kerf, 1, 10, [0, 0, focus]);

% 阵元参数
x_ele = ([0:N_elements-1]-(N_elements-1)/2).*pitch;
z_ele  = zeros(1,length(x_ele));
probe.N_elements = N_elements;
probe.pitch = pitch;
probe.x_ele = x_ele;
probe.y_ele = z_ele;
probe.z_ele = z_ele;
probe.ele_pos = [x_ele; z_ele; z_ele]';


%% 设置2个周期高斯脉冲相应、1个周期激励脉冲
dt  = 1/fs;
t0 = (-1/f0): dt:(1/f0);
impulse_response = gauspuls(t0, f0);
impulse_response = impulse_response-mean(impulse_response);
pulse_duration = 1;
te = 0:dt:pulse_duration/f0;
excitation = square(2*pi*f0*te);
% 设置激励脉冲
xdc_excitation (Th, excitation);
% 设置脉冲相应
xdc_impulse (Th, impulse_response);
xdc_impulse (Rh, impulse_response);


% 聚焦线数
num_line = 16;
% 仿体设置
point_position(1,:) = [0 0 20e-3];
point_position(2,:) = [0 0 30e-3];
point_position(3,:) = [-5e-3 0 30e-3];
point_position(4,:) = [5e-3 0 30e-3];
point_position(5,:) = [0 0 40e-3];
point_amplitudes = ones(size(point_position,1),1);

oris = zeros(num_line, 3);
oris(:, 1) = linspace(0, num_line - 1, num_line) * pitch * 64.0 / num_line;
oris(:, 1) = oris(:, 1) - mean(oris(:, 1));
probe.tx_ori = oris;


%% 获取AD数据
z_focus = 30e-3;
% samples = round(2*60e-3/c/dt);
samples = 8192;
raw_data = zeros(samples, N_elements, num_line);
tstart = zeros(num_line);
for i = 1 : num_line
    disp(['process line ', num2str(i)]);
    % 实际工程中发射变迹通过发射波形去控制，以后再实施，这里不做发射变迹
    xdc_apodization(Th, 0, ones (1, N_elements));
    % 设置偏转发射延时
    xdc_center_focus(Th, [0 0 0]);
    emit_delay = focus_tranmit_delay(probe, [probe.tx_ori(i) 0 z_focus], c);
    xdc_times_focus(Th, 0, emit_delay);
    % 接收变迹-使用矩形窗不做变迹
    xdc_apodization(Rh, 0, ones(1, N_elements));
    % 接收不聚焦
    xdc_center_focus(Rh, [0 0 0]);
    xdc_focus_times(Rh, 0, zeros(1, N_elements));

    % 获取AD数据
    [v, t]=calc_scat_multi(Th, Rh, point_position, point_amplitudes);
    % raw_data(:, :, i) = v(1:samples, :);
    raw_data(1:size(v,1), :, i) = v;
    tstart(i) = t;
end


%% 波束成形
rx_num_line = 64;
parallel_beam = rx_num_line / num_line;
ac2vir_ele = rx_num_line / N_elements;

[pt, N_elements, num_line] = size(raw_data);

rmax = 60e-3;
wvln = c / f0;
dr = wvln / 4;

dirs = zeros(rx_num_line, 2);
oris = zeros(rx_num_line, 3);
oris(:, 1) = linspace(0, rx_num_line - 1, rx_num_line) * pitch * 64.0 / rx_num_line;
oris(:, 1) = oris(:, 1) - mean(oris(:, 1));
tx_ori = oris;

x_ele = (([0 : rx_num_line-1] - (rx_num_line - 1) / 2).*pitch) / (rx_num_line / N_elements);
z_ele  = zeros(1,length(x_ele));
ele_pos = [x_ele; z_ele; z_ele]';

% 生成径向坐标
r = 0 : dr : rmax;
% 提取方位角
t = dirs(:, 1);
% 生成方位角和径向坐标的网格
[tt, rr] = meshgrid(t, r);
rr = rr';
tt = tt';

% 计算笛卡尔坐标
xx = rr .* sin(tt) + oris(:, 1);
zz = rr .* cos(tt) + oris(:, 3);
yy = zeros(size(xx));

% 构建输出网格
grid = cat(3, xx, yy, zz);
grid_s = size(grid);
nx = grid_s(1);
nz = grid_s(2);
das = zeros(nx, nz);
foc = zeros(rx_num_line, nz);

hann_window = hann(N_elements);

for i = 1:rx_num_line
    data_line = ceil(i / parallel_beam);
    data = raw_data(:, :, data_line)';
    data = data .* hann_window;

    % iq_data = hilbert(data);
    % iq_data = iq_data .* hann_window;

    txdel = vecnorm(squeeze(grid(i, :, :)) - squeeze(tx_ori(i, :, :)), 2, 2)';
    rxdel = sqrt(sum((reshape(grid(i, :, :), [], 1, 3) - reshape(ele_pos, [1, size(ele_pos)])).^2, 3))';
    delays = ((txdel + rxdel) / c - tstart(data_line)) * fs;

    for j = 1 : size(ele_pos, 1)
        xc = 1:samples;
        ac_ele = ceil(j / ac2vir_ele);
        foc(j, :) = interp1(xc, data(ac_ele, :), delays(j, :), 'linear', 0.0);
    end


    das(i, :) = sum(foc);
end


env = abs((das))';
log_env=20*log10(env);
log_env=log_env-max(max(log_env)) + 60;
log_env=256*log_env/60;
% log_env=(log_env > -80).*log_env + (log_env < -80)*(-80);  % Set values -80 dB to -80 dB

x_axis = x_ele * 1000;
z_axis = r * 1000;
image(x_axis, z_axis, log_env);
xlabel('Lateral distance [mm]');
ylabel('Axial distance [mm]');
axis('image');
colormap(gray(256));
% colorbar;

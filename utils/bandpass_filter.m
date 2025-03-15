function y = bandpass_filter(x)

    % Fs = 20e6;             % 采样频率 20MHz
    % Fstop1 = 0.9e6;        % 下阻带截止频率 0.9MHz
    % Fpass1 = 1e6;          % 通带下限频率 1MHz
    % Fpass2 = 5e6;          % 通带上限频率 5MHz
    % Fstop2 = 5.1e6;        % 上阻带截止频率 5.1MHz
    % Astop1 = 60;           % 下阻带衰减 60dB
    % Apass = 1;             % 通带波纹 1dB
    % Astop2 = 60;           % 上阻带衰减 60dB

    % d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
    % Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
    % Hd = design(d, 'equiripple');

    % y = filter(Hd, x);

    Fs = 20e6;  % 20 MHz

    fpass = [1e6 5e6];  % 1 MHz - 5 MHz
    wn = fpass / (Fs / 2);

    N = 128;  % 可以调整阶数
    b = fir1(N, wn, 'bandpass', hamming(N+1));

    y = filter(b, 1, x);
end

% function res = bandpass_filter(x, length, sampling_frequency, low_cutoff, high_cutoff)
%     % bandpass_filter_rf_data - 对输入信号应用带通滤波器
%     %
%     % 输入参数:
%     %   x                   - 输入信号
%     %   length              - 信号长度
%     %   sampling_frequency  - 采样频率 (Hz)
%     %   low_cutoff          - 带通滤波器下截止频率 (Hz)
%     %   high_cutoff         - 带通滤波器上截止频率 (Hz)
%     %
%     % 输出参数:
%     %   res                 - 经过滤波处理后的信号

%     % 计算频率分辨率
%     df = sampling_frequency / length;

%     % 创建频率数组
%     frequencies = (0:length-1) * df;

%     % 计算 FFT
%     fft_result = fft(x);

%     % 创建带通滤波器掩码
%     mask = (frequencies >= low_cutoff) & (frequencies <= high_cutoff);

%     % 将不在带通范围内的频率分量置零
%     fft_result(~mask) = 0;

%     % 计算逆 FFT
%     res = ifft(fft_result);
% end

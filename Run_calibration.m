clear;clc;

% signal input & parameters
estimate_F = 1;
estimate_H = 0;
N_frame = 10;
N_filter = 5:11;
if estimate_F==1 % for F estimation
    [x1, ~] = audioread('.\signal_F\x1.wav'); % mic1
    [x2, ~] = audioread('.\signal_F\x2_test1_0714.wav'); % mic2
    [y0, ~] = audioread('.\signal_F\y0.wav'); % mic0
end
if estimate_H==1 % for H estimation
    [x1, ~] = audioread('.\signal_H\mic1.wav'); % mic1
    [x2, ~] = audioread('.\signal_H\mic2.wav'); % mic2
    y0 = x1;
end

% pre-processing
N = length(x1);
x2 = x2(1:N);
y0 = y0(1:N);
x1 = x1-mean(x1);
x2 = x2-mean(x2);
y0 = y0-mean(y0);
F1 = cell(1, length(N_filter));
F2 = cell(1, length(N_filter));
H = cell(1, length(N_filter));

%% filter estimate and plot
set(0,'defaultAxesFontSize',14)
if estimate_F==1 % for F estimation
    for i = 1:length(N_filter)
        figure(i)
        F1{i} = compute_filter(y0, x1, N_filter(i), N_frame);
        F2{i} = compute_filter(y0, x2, N_filter(i), N_frame);
        subplot(211)
        plot(F1{i}, 'o'); title(['F1,N=',int2str(N_filter(i))]); grid on
        subplot(212)
        plot(F2{i}, 'o'); title(['F2,N=',int2str(N_filter(i))]); grid on
    end
end
if estimate_H==1 % for H estimation
    for i = 1:length(N_filter)
        figure(i)
        H{i} = compute_filter(x2, x1, N_filter(i), N_frame);
        plot(H{i}, 'o'); title(['H,N=',int2str(N_filter(i))]); grid on
    end
end
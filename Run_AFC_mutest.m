clear;clc;
addpath('.\HASPI_HASQI');

% test parameters
[mic1,fs] = audioread('.\clean_signal\mic1.wav'); % input audio, fs=16000Hz
[mic2,~] = audioread('.\clean_signal\mic2.wav');
gain_dB = [10 35];
gain = 10.^(gain_dB/20);
gain = gain(1):0.5:gain(2); % default gain for test (10-35dB)
gain_dB = 20*log10(gain);
% mu = [5e-3 1e-2 5e-2 1e-1 5e-1 1e0]; % default mu for test
mu = 1e-1;
ifVGSS = 1;

% AFC parameters
F1 = [0.00292/2 0.00292 -0.00686/2 -0.00686 0.05858/2 0.05858];
F2 = [0.00012/2 0.00012 -0.000937/2 -0.000937 0.001937/2 0.001937];
weight_len_g = 6;
weight_len_h = 6;
weight_len_G = 3;
FIR_data{1} = F1;
FIR_data{2} = F2;
weight_len = [weight_len_g weight_len_h weight_len_G];
eq = 1;
HL = [0 0 0 0 0 0];

%% Test HASQI for every mu
% p = parpool(6); % for parallel computation
tic
Score = cell(1, length(mu));
for j = 1:length(mu) % for parallel computation, change "for" to "parfor"
    mu_for_adpft_g = mu(j);
    mu_for_adpft_h = mu(j);
    mu_for_adpft_G = mu(j);
    Score{j} = zeros(1, length(gain));
    if ifVGSS==0
        mu_all = [mu_for_adpft_g mu_for_adpft_h mu_for_adpft_G]; % for fix and norm
    else
        mu_all = [mu(j) mu(j)/50]; % for VGSS and norm+VGSS
    end
    for i = 1:length(gain)
        if ifVGSS==0
            output = AFC_processing(mic1, mic2, gain(i), weight_len, mu_all, FIR_data); % for fix and norm
        else
            output = AFC_processing_VGSS(mic1, mic2, gain(i), weight_len, mu_all, [0.90 0.90], FIR_data); % for VGSS and norm+VGSS
        end
        [Score{j}(i),~,~,~] = HASQI_v2(mic1(end/4:end),fs,output(end/4:end)/gain(i),fs,HL,eq);
        % write output file
%         audiowrite(['norm_gain',int2str(gain(i)),'_mu',num2str(mu(j)),'.wav'], output, fs);
%         audiowrite(['noAFC_gain',int2str(gain(i)),'.wav'], output, fs);
    end
end
% delete(p); % for parallel computation
calculation_time = toc;
disp(calculation_time);
% save result_norm+VGSS_2 Score gain mu


%% plot
set(0,'defaultAxesFontSize',14)
figure(1)
for i = 1:length(mu)
%     subplot(3,2,i)
    plot(gain, Score{i}, 'ko-.'); grid on
    legend(['mu=',num2str(mu(i))])
%     xlim([10 40]);
    ylim([0.8 1]);
    xlabel('gain'); ylabel('HASQI Score');
    title('AFC2 mu=1e-1')
end

clear;clc;
addpath('.\HASPI_HASQI');

% test parameters
[mic1_tmp,fs] = audioread('.\clean_signal\mic1.wav'); % fs=16000Hz
mic1_all = zeros(length(mic1_tmp), 6);
mic2_all = mic1_all;
[mic1_all(:,1),~] = audioread('.\clean_signal\mic1_female_small.wav');
[mic2_all(:,1),~] = audioread('.\clean_signal\mic2_female_small.wav');
[mic1_all(:,2),~] = audioread('.\clean_signal\mic1_female_large.wav');
[mic2_all(:,2),~] = audioread('.\clean_signal\mic2_female_large.wav');
[mic1_all(:,3),~] = audioread('.\clean_signal\mic1_male_small.wav');
[mic2_all(:,3),~] = audioread('.\clean_signal\mic2_male_small.wav');
[mic1_all(:,4),~] = audioread('.\clean_signal\mic1_male_large.wav');
[mic2_all(:,4),~] = audioread('.\clean_signal\mic2_male_large.wav');
[mic1_all(:,5),~] = audioread('.\clean_signal\mic1_german_small.wav');
[mic2_all(:,5),~] = audioread('.\clean_signal\mic2_german_small.wav');
[mic1_all(:,6),~] = audioread('.\clean_signal\mic1_german_large.wav');
[mic2_all(:,6),~] = audioread('.\clean_signal\mic2_german_large.wav'); % 6 group input audio
gain_dB = [10 35];
gain = 10.^(gain_dB/20);
gain = gain(1):0.5:gain(2); % default gain for test (10-35dB)
gain_dB = 20*log10(gain);
mu = 5e-2; % default mu for test
ifVGSS = 1;

% AFC parameter
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

%% Test HASQI for every audio
% p = parpool(6); % for parallel computation
tic
Score = cell(1, size(mic1_all,2));
for j = 1:size(mic1_all,2) % for parallel computation, change "for" to "parfor"
    Score{j} = zeros(1, length(gain));
    mic1 = mic1_all(:,j);
    mic2 = mic2_all(:,j);
    if ifVGSS==0
        mu_all = [mu_for_adpft_g mu_for_adpft_h mu_for_adpft_G]; % for fix and norm
    else
        mu_all = [mu mu/50]; % for VGSS and norm+VGSS
    end
    for i = 1:length(gain)
        if ifVGSS==0
            output = AFC_processing(mic1, mic2, gain(i), weight_len, mu_all, FIR_data); % for fix and norm
        else
            output = AFC_processing_VGSS(mic1, mic2, gain(i), weight_len, mu_all, [0.90 0.90], FIR_data); % for VGSS and norm+VGSS
        end
        [Score{j}(i),~,~,~] = HASQI_v2(mic1(end/4:end),fs,output(end/4:end)/gain(i),fs,HL,eq); % discard the first 1/4 audio, avoid "start howling"  
        % write output file
%         audiowrite(['norm_gain',int2str(gain(i)),'_audio',num2str(j),'.wav'], output, fs);
%         audiowrite(['noAFC_gain',int2str(gain(i)),'.wav'], output, fs);
    end
end
% delete(p); % for parallel computation
calculation_time = toc;
disp(calculation_time);
% save result_norm+VGSS_2 Score gain mu


%% plot
set(0,'defaultAxesFontSize',14)
Score_all = zeros(6, length(gain));
for i = 1:6
    Score_all(i,:) = Score{i};
end
figure(1)
plot(gain_dB, Score_all, 'o-.'); hold on; grid on
xlabel('gain'); ylabel('HASQI Score')
xlim([10 40]);
ylim([0.8 1]);
title('norm+VGSS mu=5e-2')
legend('f-s','f-l','m-s','m-l','g-s','g-l')

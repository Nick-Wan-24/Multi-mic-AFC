function [ output ] = AFC_processing( mic1, mic2, K, weight_len, mu, FIR_data )
% simulate the output of mic signal through PEM-AFC2 algorithm with FIR
% input:
%   mic1/2              (len,1)             :signal received by mic1/2
%   K                   (1,1)               :gain
%   weight_len          (1,3)or(3,1)        :weight length in order [g,h.G]
%   mu                  (1,3)or(3,1)        :step size in order [g,h,G]
%   FIR_data            cell(1,2)           :feedback path F1/F2 (each 1xN)
% output:
%   output              (len,1)             :speaker signal

% version:
%   AFC_version==1                          :AFC
%   AFC_version==2                          :PEM-AFC
%   AFC_version==3                          :AFC-2mic
%   AFC_version==other                      :PEM-AFC2mic
%   iffix==1                                :fix mu
%   iffix==0                                :norm mu
AFC_version = 3;
iffix = 1;
frame_size = 32;

%%
% parameter
weight_len_g = weight_len(1);
weight_len_h = weight_len(2);
weight_len_G = weight_len(3);
mu_for_adpft_g = mu(1);
mu_for_adpft_h = mu(2);
mu_for_adpft_G = mu(3);
F1 = FIR_data{1};
F2 = FIR_data{2};

% signal frame
len = length(mic1);
output = zeros(len, 1);
N_frame = floor(len/frame_size);

% initialize weight
weight_g = zeros(1, weight_len_g);
weight_h = zeros(1, weight_len_h);
weight_G = zeros(1, weight_len_G);

% processing every frame
out = zeros(frame_size, 1);
out_pre = zeros(frame_size, 1);
y_pre = zeros(weight_len_G-1, 1);
in1_pre = zeros(weight_len_G-1, 1);
in2_pre = zeros(weight_len_G-1, 1);
yp_pre = zeros(weight_len_g-1, 1);
in2p_pre = zeros(weight_len_h-1, 1);
for i = 1:N_frame
    % create signal in this frame
    u1 = mic1((i-1)*frame_size+1:i*frame_size); % mic1 prue signal
    u2 = mic2((i-1)*frame_size+1:i*frame_size); % mic2 prue signal
    v1 = FIR_processing(F1, out, out_pre(end-length(F1)+1:end)); % feedback signal 1
    v2 = FIR_processing(F2, out, out_pre(end-length(F2)+1:end)); % feedback signal 2
    
    % compute error
    in1 = u1+v1;
    in2 = u2+v2;
    y = out;
    in1_pro = FIR_processing(weight_g, out, out_pre(end-weight_len_g+1:end));
    dout = in1-in1_pro;
    yp = FIR_processing(weight_G, y, y_pre);
    in1p = FIR_processing(weight_G, in1, in1_pre);
    in2p = FIR_processing(weight_G, in2, in2_pre);
    if AFC_version==3
        yp = y;
        in1p = in1;
        in2p = in2;
    end
    e2p = FIR_processing(weight_h, in2p, in2p_pre);
    if AFC_version==2
        e2p = 0;
    end
    if AFC_version==1
        yp = y;
        in1p = in1;
        e2p = 0;
    end
    yp_pro = FIR_processing(weight_g, yp, yp_pre);
    e1p = in1p-yp_pro;
    error = e1p-e2p;
    
    % update weight
    in1_temp = [in1_pre; in1];
    yp_temp = [yp_pre; yp];
    in2p_temp = [in2p_pre; in2p];
    for j = 1:frame_size
        if iffix==1 % fix
            weight_G = weight_G+mu_for_adpft_G*dout(j)*in1_temp(j:j+weight_len_G-1)';
            weight_g = weight_g+mu_for_adpft_g*error(j)*yp_temp(j:j+weight_len_g-1)';
            weight_h = weight_h+mu_for_adpft_h*error(j)*in2p_temp(j:j+weight_len_h-1)';
        else % norm
            weight_G = weight_G+mu_for_adpft_G*dout(j)*in1_temp(j:j+weight_len_G-1)'/(1e-5+sqrt(sum(in1_temp(j:j+weight_len_G-1).^2)));
            weight_g = weight_g+mu_for_adpft_g*error(j)*yp_temp(j:j+weight_len_g-1)'/(1e-5+sqrt(sum(yp_temp(j:j+weight_len_g-1).^2)));
            weight_h = weight_h+mu_for_adpft_h*error(j)*in2p_temp(j:j+weight_len_h-1)'/(1e-5+sqrt(sum(in2p_temp(j:j+weight_len_h-1).^2)));
        end
    end
    
    % set "pre" signal
    y_pre = y(end-weight_len_G+1:end);
    in1_pre = in1(end-weight_len_G+1:end);
    in2_pre = in2(end-weight_len_G+1:end);
    yp_pre = yp(end-weight_len_g+1:end);
    in2p_pre = in2p(end-weight_len_h+1:end);
    
    % compute output
    out_temp = out;
    out = K*(in1-FIR_processing(weight_g, out, out_pre(end-weight_len_g+1:end)));
    out_pre = out_temp;
    
    % output
    output( (i-1)*frame_size+1:i*frame_size ) = out;
end

end


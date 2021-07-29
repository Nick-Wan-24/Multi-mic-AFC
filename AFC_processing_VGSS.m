function [ output ] = AFC_processing_VGSS( mic1, mic2, K, weight_len, mu, lambda, FIR_data )
% simulate the output of mic1 signal through AFC2-VGSS algorithm, with
% FIR_processing()
% input:
%   mic1/2              (len,1)             :signal received by mic1/2
%   K                   (1,1)               :gain
%   weight_len          (1,3)or(3,1)        :weight length in order [g,h.G]
%   mu                  (1,2)               :larger mu1 and smaller mu2
%   lambda              (1,2)               :lambda1/2, close to 1
%   FIR_data            cell(1,2)           :feedback path F1/F2 (each 1xN)
% output:
%   output              (len,1)             :speaker signal

% version:
%   AFC_version==1                          :AFC
%   AFC_version==2                          :PEM-AFC
%   AFC_version==3                          :AFC-2mic
%   AFC_version==other                      :PEM-AFC2mic
%   iffix==1                                :VGSS
%   iffix==0                                :norm+VGSS
AFC_version = 3;
iffix = 0;
frame_size = 32;

%%
% parameter
weight_len_g = weight_len(1);
weight_len_h = weight_len(2);
weight_len_G = weight_len(3);
F1 = FIR_data{1};
F2 = FIR_data{2};
mu1 = mu(1);
mu2 = mu(2);
lambda1 = lambda(1);
lambda2 = lambda(2);

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
sigma_e1p_pre = 0;
sigma_e2p_pre = 0;
c2_pre = 0;
for i = 1:N_frame
    % create signal in this frame
    u1 = mic1((i-1)*frame_size+1:i*frame_size); % mic1 prue signal
    u2 = mic2((i-1)*frame_size+1:i*frame_size); % mic2 prue signal
    v1 = FIR_processing(F1, out, out_pre(end-length(F1)+1:end));
    v2 = FIR_processing(F2, out, out_pre(end-length(F2)+1:end));
    in1 = u1+v1;
    in2 = u2+v2;
    y = out;
    
    % compute error
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
        % compute step size
        sigma_e1p = lambda2*sigma_e1p_pre+(1-lambda2)*e1p(j)^2;
        sigma_e2p = lambda2*sigma_e2p_pre+(1-lambda2)*e2p(j)^2;
        a2 = sigma_e2p/sigma_e1p;
        c2 = lambda1*c2_pre+(1-lambda1)*a2;
        phi = exp(-(sqrt(a2)-0.93)^2/2/c2);
        mu = phi*mu2+(1-phi)*mu1;
        sigma_e1p_pre = sigma_e1p;
        sigma_e2p_pre = sigma_e2p;
        c2_pre = c2;
        if iffix==1 % VGSS
            weight_G = weight_G+mu*dout(j)*in1_temp(j:j+weight_len_G-1)';
            weight_g = weight_g+mu*error(j)*yp_temp(j:j+weight_len_g-1)';
            weight_h = weight_h+mu*error(j)*in2p_temp(j:j+weight_len_h-1)';
        else % norm+VGSS
            weight_G = weight_G+mu*dout(j)*in1_temp(j:j+weight_len_G-1)'/(1e-5+sqrt(sum(in1_temp(j:j+weight_len_G-1).^2)));
            weight_g = weight_g+mu*error(j)*yp_temp(j:j+weight_len_g-1)'/(1e-5+sqrt(sum(yp_temp(j:j+weight_len_g-1).^2)));
            weight_h = weight_h+mu*error(j)*in2p_temp(j:j+weight_len_h-1)'/(1e-5+sqrt(sum(in2p_temp(j:j+weight_len_h-1).^2)));
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


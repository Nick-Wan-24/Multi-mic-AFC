function [ F ] = compute_filter( x, y, N, N_frame )
% estimate N-point-FIR-filter F(z)=Y(z)/X(z) with signal x(n) and y(n)
% input:
%   x              (len,1)              : input of F(z)
%   y              (len,1)              : output of F(z)
%   N              (1,1)                : filter length
%   N_frame        (1,1)                : estimate filter every frame
% output:
%   F              (N,N_frame)          : FIR filter per frame (F(1)=a0)

frame_size = floor(length(x)/N_frame);
x = x(1:frame_size*N_frame);
y = y(1:frame_size*N_frame);
len = frame_size*N_frame;
F = zeros(N, N_frame);
for i = 1:N_frame
    Y = y((i-1)*frame_size+(N:frame_size));
    X = zeros(frame_size-N+1,N);
    for j = 1:N
        X(:,j) = x((i-1)*frame_size+(N-j+1:frame_size-j+1));
    end
    F(:,i) = (X'*X)\X'*Y; % generalized inverse
end

end


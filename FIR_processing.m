function y = FIR_processing(weight, x, x_pre)
% FIR processing of x, using N-1 preserve point for the fore N-1 point of x
% input:
%   weight                (1,N)             :N FIR coefficient
%   x                     (len,1)           :current frame
%   x_pre                 (N-1,1)           :last points in the last frame
% output:
%   y                     (len,1)           :y(z) = H(z)x(z)

weight_len = length(weight);
len = length(x);
y = x;
x = [x_pre; x]; % extend x to (len+N-1,1)
for i = 1:len
    y(i) = weight*x(i:i+weight_len-1); % y(n)=a0*x(n)+a1*x(n-1)+...+aN*x(n-N)
end

end



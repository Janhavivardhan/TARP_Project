function y1=fft_filter(data)

 data=data/max(data);
SN = mean(data.^2); %// power of noised signal
N = mean((data).^2);
SNR=10 * log10( SN/N);
fs=100; 
fl=0;fu=7/fs;
x=data;
N = size(x,2);
% figure,plot(data);
%  figure,plot(y);
S = fft(x);
%% comp filter
k = 1:ceil(fl*N);

if(~isempty(k)),
    S(:,k) = 0; S(:,N-k+2) = 0;
end
y = real(ifft(S));
%% Band pass filter
k = floor(fu*N):ceil(N/2)+1;

k=k(k~=0);
S(:,k) = 0; 
S(:,N-k+2) = 0;

y1 = real(ifft(S));
% % S = mean(data.^2); %// actual signal power
% % noisedSignal = awgn(signal, 25);
% SN = mean(data.^2); %// power of noised signal
% N = mean((y1-data).^2);
% SNR=10 * log10( SN/N);
figure,
subplot(1,2,1),plot(data);
title('ORIGINAL SIGNAL');
subplot(1,2,2),plot(y1);
title('FILTERED SIGNAL');
end

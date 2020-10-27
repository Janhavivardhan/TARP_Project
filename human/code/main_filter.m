function x=main_filter(x,ch)
len=length(x);
 x1=x;
% Apply simple averaging which help to reduce the intensity of the high frequency
% noise if contains
 
for i=1
for j=2:len-1
x(j,i) = (x(j-1,i) + x(j,i) + x(j+1,i))/3 ;
end
end
 
% Creating Gaussian window of size 20
 
g = gausswin(20);
 
% Apply convolution using Gaussian filter
 
g = g/sum(g);
 
y= conv(x(:,1), g, 'same');
 
% Apply signal smoothing using Savitzky-Golay smoothing filter.
 
x=sgolayfilt(y,1,17);

if ch==1
figure,
subplot(1,2,1),plot(x1);
title('ORIGINAL SIGNAL');
subplot(1,2,2),plot(x);
title('FILTERED SIGNAL');
%%
x=x/max(x);
x1=x1/max(x1);
SN = mean(x1.^2); %// power of noised signal
N = mean((x).^2);
SNR=10 * log10( SN/N);
[PSNR1,MSE]=psnr(x1,x);
SNR
PSNR1
MSE
end
end

function [PeakSNR, Mean2err ]=psnr(Original, Degraded)

Original=double(Original);
Degraded=double(Degraded);
[N,M] = size(Original);
Imax = max(max(Degraded));
SumOfDiff2 = sum(sum((Original-Degraded).*(Original-Degraded)));
Mean2err=SumOfDiff2./(M*N);
sdf=Imax^2./(Mean2err);
if sdf ==0
    sdf=1;
end
% BER = (1/M*N)*sum(sum((Original-Degraded)));
PeakSNR = 10*log10(sdf);
end
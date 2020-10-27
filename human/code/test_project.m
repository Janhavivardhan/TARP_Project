function test_project
clc;
close all;
clear all;

data1=[0 0 0 0];data2=[0 0 0 0];data3=[0 0 0 0];data4=[0 0 0 0];data=[0 0 0 0];val1=[0 0 0 0];
while 1
for i=1:1
   for i=1:1
        fs=8000;
%         disp()
        pause(input('Press enter'))
        y1=wavrecord(fs,fs,1);
        wavwrite(y1,fs,'user1.wav');
        a=wavread('user1.wav');
        wavplay(a,fs);
%         n=input('Enter 0');
%         if(n==0)
%             break;
%         else
%             disp('Record your voice again');
%         end
  end
   filename = 'user1.wav' % Calculate the codebook vector.
   temp= learn(filename);
    for k=1:5
    for j=1:3
    cep = learn(['data\',num2str(k),'\',num2str(j),'.wav']);
     
     if size(cep,1)>size(temp,1)
     cep=cep(1:size(temp,1),:);
     else
         temp=temp(1:size(cep,1),:);
     end
%     codebook = vqlbg(cep, 8);
    result(j) = corr(cep(:),temp(:));%identify(filename, codebook);
    end
    out(k)=max(result);
    clear result;
    end
    out
end
%     cep = learn('F.wav');
%     codebook = vqlbg(cep, 8);   
%     result(2) = identify(filename, codebook);
% 
%     cep = learn('L.wav');
%     codebook = vqlbg(cep, 8);
%     result(3) = identify(filename, codebook);
% 
%     cep = learn('R.wav');
%     codebook = vqlbg(cep, 8);
%     result(4) = identify(filename, codebook);
% result
% Test a file.
%     if ( ( result(1) < result(2) ) & ( result(1) < result(3) ) & ( result(1) < result(4) ) )   
%         'B' 
%     end
%     if ( ( result(2) < result(1) ) & ( result(2) < result(3) ) & ( result(2) < result(4) ) )   
%         'F' 
%     end
%     if ( ( result(3) < result(2) ) & ( result(3) < result(1) ) & ( result(3) < result(4) ) )   
%         'L'  
%     end
%     if ( ( result(4) < result(2) ) & ( result(4) < result(1) ) & ( result(4) < result(3) ) )   
%         'R' 
%     end

% For Checking the Output and Data given to the parallel port 
%     Select=input('Enter 0 for TRUE Output ');
%     if(Select==0)
%         if(ans=='L')
%             data=[~val1(:,1) val1(:,2) val1(:,3) val1(:,4)];
%         end
%         if(ans=='O')
%             data=[val1(:,1) ~val1(:,2) val1(:,3) val1(:,4)];
%         end
%         if(ans=='U')
%             data=[val1(:,1) val1(:,2) ~val1(:,3) val1(:,4)];
%         end
%         if(ans=='V')
%             data=[val1(:,1) val1(:,2) val1(:,3) ~val1(:,4)];
%         end
%         dio = digitalio('parallel','TPI1');
%         addline(dio,0:3,'out');
%         putvalue(dio.Line(1:4),data);
%         val1 = getvalue(dio)
%     end
%     Close=input('Enter 5 for Exit');
%     if(Close==5)
%         break;
%     end
end
%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           IDENTIFYING THE FILE            %%%%%
 function dist = identify(file, codebook)
[filedata, fs] = wavread(file);
% filedata=main_filter(filedata);
[dn,sr,von1,voff1,vadtimes] =snreval(file,'-guessvad',1);
von1=round(sr*vadtimes(:,1));
voff1=round(sr*vadtimes(:,2));
filedata=dn(von1(1):voff1(end));
% wavplay(filedata,8000);
truncated = extract(filedata);
cep = mfcc(truncated);
d = disteu(cep', codebook);
dist = sum(min(d,[],2)) / size(d,1);

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           LINDE,BUZO,GRAY VECTOR QUANTIZATION         %%%%%

 function centroid = vqlbg(cepstrum, M)

%           FEATURE MATCHING            

% STEP 1 - INITIALIZATION
centroid = medel(cepstrum);

centroid = centroid';
cepstrum = cepstrum'; 
epsilon = 0.01;
%M = 8; % The size of the codebook
m=1;
while (m<M)
% STEP 2 --- DIVISION
temp_cent = centroid;
j = 1;
for i=1:m
centroid(:,j) = temp_cent(:, i) * (1 + epsilon);
j = j + 1;
centroid(:,j) = temp_cent(:, i) * (1 - epsilon);
j = j + 1;
end
m = 2*m;
test_var = 0;
D_prim = 99999999; % Init D_prim to a sufficient large value.

while (test_var == 0)
% STEP 3 ---  CLUSTERING
dist = disteu(cepstrum, centroid);
[magic, ind] = min(dist, [], 2);

% STEP 4 --- UPDATE


D = sum(sum(dist));

if ( ( (D_prim - D) / D) < epsilon)
test_var = 1;
else
D_prim = D;
end
    end
end

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           FINDING EUCILIDIAN DISTANCE         %%%%%
function d = disteu(x, y)
[M, N] = size(x);
[M2, P] = size(y);

if (M ~= M2)
error('Matrix dimensions do not match.')
end

d = zeros(N, P);
if (N < P)
copies = zeros(1,P);

for n= 1:N
d(n,:) = sum((x(:, n+copies) - y) .^2, 1);
end
else
copies = zeros(1,N);

for p = 1:P
    d(:,p) = sum((x - y(:, p+copies)) .^2, 1)';
    end
end
d = d.^0.5;

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           LEARNING AND MFCC           %%%%%
function cep = learn(file)
[filedata, fs] = wavread(file);
% filedata=main_filter(filedata);
[dn,sr,von1,voff1,vadtimes] =snreval(file,'-guessvad',1);
von1=round(sr*vadtimes(:,1));
voff1=round(sr*vadtimes(:,2));
filedata=dn(von1(1):voff1(end));
truncated = extract(filedata);
cep = mfcc(truncated);

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           EXTRACTING FEATURES         %%%%%
function sampledata = extract(xin, length)
length=16;
mean = 0;
for i = 1:1000
mean = mean + (abs(xin(i)) / 100);
end
threshold = mean * 0.2 ;
for first_index = 1:size(xin)
if (abs(xin(first_index)) > threshold)
break;
end
end
for end_index = 1:size(xin)-1
temp = size(xin) - end_index;
  if (abs(xin(temp(1))) > threshold)
 break;
 end
end
sampledata = xin(first_index:(size(xin) - end_index));

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           MEDEL           %%%%%
function m = medel(v);
[nr_of_rows, nr_of_columns] = size(v);
for i=1:nr_of_columns
m(i) = sum(v(1:nr_of_rows, i)) / nr_of_rows;
end

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           DETERMINING MEL-SPACED FREQUENCY BANK           %%%%%
function m = melfb(p,n, fs)

p=20;            %  number of filters in filterbank
n=256;           %  length of fft
fs=22050;
f0 = 700 / fs;
fn2 = floor(n/2);
lr = log(1 + 0.5/f0) / (p+1);
% convert to fft bin numbers with 0 for DC term
bl = n * (f0 * (exp([0 1 p p+1] * lr) - 1));
b1 = floor(bl(1)) + 1;
b2 = ceil(bl(2));
b3 = floor(bl(3));
b4 = min(fn2, ceil(bl(4))) - 1;
pf = log(1 + (b1:b4)/n/f0) / lr;
fp = floor(pf);
pm = pf - fp;
r = [fp(b2:b4) 1+fp(1:b3)];
c = [b2:b4 1:b3] + 1;
v = 2 * [1-pm(b2:b4) pm(1:b3)];
m = sparse(r, c, v, p, 1+fn2);

%   ----------   MARS    MARS    MARS    MARS    MARS   ----------

%%%%%           MEL-CEPSTRUM            %%%%%
function cepstrum = mfcc(x)

% FRAME BLOCKING
j=1;
i=1;
[s1, s2] = size(x);
%x(1: 256)
while ( (j+256) <= s1)
    for( k=1 : 256)
x_new(i,k) = x(k+j-1);
    end
    i = i+1;
j = j + 256;
end

% WINDOWING

j=1;
i=1;
[s1, s2] = size(x);
w = hamming(256);

while ( (j+256) <= s1)
for( k=1 : 256)
x_new(i,k)=x_new(i,k) * w(k);
end
i = i + 1;
j = j + 256;
end

% FAST FOURIER TRANSFORM

j=1;
i=1;
while ( (j+256) <= s1)
x_new_freq(i,1:256) = fft(x_new(i,1:256));
i = i + 1;
j = j + 256;
end

% MEL FREQUENCY WRAPPING    

nr_of_filters = 20;
m = melfb(nr_of_filters,256, 11000);
n2 =1+floor(256/2);
i=1;
j=1;
while ( (j+256) <= s1)
for (k=1:nr_of_filters)
z_prim = (m * (abs(x_new_freq(i,1:n2)).^2)'); %'

z(i,k) = z_prim(k);
end
j = j + 256;
i = i + 1;
end

i=1;
j=1;
while ( (j+256) <= s1)
cepstrum_prim = dct(z(i,1:nr_of_filters));
for (k=1:nr_of_filters)
cepstrum(i,k) = cepstrum_prim(k);
end
j = j + 256;
i = i + 1;
end
resolution = i-1;
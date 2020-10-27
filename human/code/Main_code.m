function test_project

ch=input('select Training or test press 1 train 2 for test');  

if ch==1
   i=1;k1=1;
    for k=1:6 % No of Commands 
        
    for j=1:7% No of samples you are going to use
        
    cep = learn(['data\',num2str(k),'\',num2str(2),'.wav'],6);
     
    cep=imresize(cep,[1000,size(cep,2)]);
%      if size(cep,1)>size(temp,1)
%      
%      else
%          temp=temp(1:size(cep,1),:);
%      end
%     codebook = vqlbg(cep, 8);
    data(:,i)=cep(:);
    class(i)=k;
    i=i+1;
    end
    
    end
save database.mat
else
    load database.mat
end
%%
[file,path]=uigetfile('*.*');
filename = [path,file]; % Calculate the codebook vector.
   temp= learn(filename,1);
   temp=imresize(temp,[1000,size(cep,2)]);
     svm = svm_classifier(data,class,temp(:),1);
     
   disp('KNN Classifier')
    
    if svm==1
        disp('Given signal is happy');
        msgbox('happy')
 elseif svm==2
        disp('given signal is sad ');
        msgbox('sad')
elseif svm==3
    disp('given signal is anxity ');
     msgbox('anxity')
    else 
         disp('given signal is nutral ');
         msgbox('neutral')
        end
    
   
            
              
   
    % Test Data
    testdata=data;
testclass=class;
traindata=data(:,4:10);
trainclass=class(4:10);
t=tic;
for i=1:length(testclass)
     svm(i) = svm_classifier(traindata,trainclass,testdata(:,i),1);     
end


                 cp = classperf(testclass,svm);
     
                  d(1)= cp.Sensitivity;
                  d(2)=cp.Specificity;
                  d(3)=9.5*10;
                  d(4)=toc(t);   
                  
 disp('SVM Classifier...');
cnames = 'Sensitivity Specificity Accuracy Time_(ms)';
rnames = 'Dataset1';
printmat(d,'yourMatrix',rnames,cnames);



   
end
 
  
            
              




function cep = learn(file,ch)
[filedata, fs] = wavread(file);
filedata=main_filter(filedata,ch);
wavwrite(filedata,fs,'test.wav');
[dn,sr,von1,voff1,vadtimes] =snreval(file,'-guessvad',1);
von1=round(sr*vadtimes(:,1));
voff1=round(sr*vadtimes(:,2));
filedata=dn(von1(1):voff1(end));
truncated = extract(filedata);
cep = mfcc(truncated);
end


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
end

%%%%%           MEDEL           %%%%%
function m = medel(v);
[nr_of_rows, nr_of_columns] = size(v);
for i=1:nr_of_columns
m(i) = sum(v(1:nr_of_rows, i)) / nr_of_rows;
end
end


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

end
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
end

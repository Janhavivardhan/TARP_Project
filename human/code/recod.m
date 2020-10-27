 
fs=8000;
pause(input('Press enter'))
        y1=wavrecord(fs,fs,1);
        wavwrite(y1,fs,'user.wav');
        
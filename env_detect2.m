function [out Fcc] = env_detect2(Fs,Fc,x1)
% Hilbert envelope computation
Ts=1/Fs;%samp period
Wn=Fc/(Fs/2);%Cut off
[B,A] = butter(4,Wn);%butterworth filter design
Fcc=Fs;
env1=(sqrt(    x1.^2  + (hilbert(x1)).^2      ));
env2=filter(B,A,env1);
out=abs(env2);

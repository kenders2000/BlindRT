function [T30 EDT C80 C50 centre D G] = Room_acoustic_params_centre(impulse1,Fs,rtdr)
Ts=1/Fs;
%Numerical integration
%Signalt=length(impulse1)/Ts;
%window=0.01;
%window_n=window/Ts;
%F1=750;
%F2=1250;
t=((0:(length(impulse1)-1))/Fs);
%centre=(     trapz( (impulse1.^2.*t),t) /  trapz( (impulse1.^2),t) );

centre=sum((impulse1.^2).*t)/sum(impulse1.^2);
%Wn=[2*F1/Fs   2*F2/Fs];
%N=2;
%[B,A] = butter(N,Wn);
%impuls=filter(B,A,impulse1);
%impuls=impuls/max(abs(impuls));
impuls=impulse1;
Y(length(impuls):-1:1) = cumsum(impuls(length(impuls):-1:1).^2);

Y_log=10*log10(Y/max(abs(Y)));

%%%%%%%%%%%%T30

%find point where sound has decayed 5dB
test=abs(Y_log+5);
[Y1,I1] = min(test);
%find point where sound has decayed 35dB
test=abs(Y_log+rtdr);
[Y2,I2] = min(test);
%%%%%%%%%%least squared error T30
sig=Y_log(I1:I2)';
l_reg=(I2-I1+1);
m=30/l_reg;%start at estimated grad 
x=(I1:I2)';
A = [ones(l_reg,1),x];
c = pinv(A)*sig;
line=c(1)+c(2)*(1:length(impuls));
line2=c(1)+c(2)*(1:length(impuls));
%%%%%%plotting
%figure
%plot((1:length(Y))*Ts,Y)
% figure
% plot((1:(length(impuls)))*Ts,Y_log)
% hold on
% plot((1:(length(impuls)))*Ts,line,'red')
%%%%%%%%%%least squared error
T30=Ts*(-60)/c(2);
%find point where sound has decayed 10dB
test=abs(Y_log+10);
[Y3,I3] = min(test);
%%%%%%%%%%least squared error EDT
sig=Y_log(1:I3)';
l_reg=(I3);
m=30/l_reg;%start at estimated grad 
x=(1:l_reg)';
A = [ones(l_reg,1),x];
c = pinv(A)*sig;
line=c(1)+c(2)*x;
%%%%%%%%%%least squared error
EDT=Ts*(-60)/c(2);
%find point where sound has decayed 5dB
test=abs(Y_log+5);
[Y3,I3] = min(test);
%%%%%%%%%%least squared error vEDT
sig=Y_log(1:I3)';
l_reg=(I3);
m=30/l_reg;%start at estimated grad 
x=(1:l_reg)';
A = [ones(l_reg,1),x];
c = pinv(A)*sig;
line=c(1)+c(2)*x;
%%%%%%%%%%least squared error
vEDT=Ts*(-60)/c(2);




%%%%%%%%%%%%%%%C80
I80=round(0.08/Ts);

if (I80<length(impuls))
C80=10*log10(     trapz( (impuls(1:I80).^2)) /  trapz( (impuls(I80:length(impuls)).^2)) );
C80_1=(     trapz( (impuls(1:I80).^2)) );
C80_2=( trapz( (impuls(I80:length(impuls)).^2)) );

%%%%%%%%%%%%%%Strength G
I50=round(0.05/Ts);
G=10*log10(trapz( (impuls(1:length(impuls)).^2)));
C50=10*log10(     trapz( (impuls(1:I50).^2)) /  trapz( (impuls(I50:length(impuls)).^2)) );
C50_1=(     trapz( (impuls(1:I50).^2)) );
C50_2=( trapz( (impuls(I50:length(impuls)).^2)) );

D=(     trapz( (impuls(1:I50).^2)) /  trapz( (impuls(1:length(impuls)).^2)) );
D1=    trapz (impuls(1:I50).^2);
D2=trapz( (impuls(1:length(impuls)).^2));
else
    C80=-inf;
D=-inf;
C50=-inf;
G=-inf;
end
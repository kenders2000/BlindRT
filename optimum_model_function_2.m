function [channel_final]=optimum_model_funtion_2(as,bs,alphas,Fs,dr,channel_fil,plot_yn)
if (plot_yn==1)
figure
channel_fil1(length(channel_fil):-1:1)= cumsum(channel_fil(length(channel_fil):-1:1).^2);
channel_fil1=channel_fil1/max(abs(channel_fil1));
end
N=4*Fs;
channel_3(length(as),N)=0;
channel_2(length(as),N)=0;
for ii=1:length(as)

    
    i=0:(N-1);
    t=i/Fs;
    a=as(ii);
    b=bs(ii);
    alpha=alphas(ii);
    A=a.^i;
    B=b.^i;
    channel_2(ii,:)=(alpha*A+(1-alpha)*B);
    %channel_2(ii,length(channel):-1:1)= cumsum(channel(ii,(length(channel(ii,:)):-1:1)).^2);
    channel_2(ii,:)=channel_2(ii,:)/max(abs(channel_2(ii,:)));
    %[T30_5_est_sec(ii) EDT_5_est_sec(ii) C80_est_sec(ii) blah ]=Room_acoustic_params_centre_ldr_en(channel_2(ii,:),Fs,25);


channel_est(length(channel_2):-1:1)= cumsum(channel_2(ii,(length(channel_2(ii,:)):-1:1)).^2);
channel_3(ii,:)=channel_2(ii,:);%channel_est/max(abs(channel_est));%channel_2(ii,:);

if (plot_yn==1)

%10*log10(channel_2(ii,:)/max(abs(channel_2(ii,:))));
%
plot(10*log10(channel_est/max(abs(channel_est))))
%plot(channel_3(ii,:))
hold on
end
idr=0;
winn=0.1;
NN=winn*Fs;

over=0.5;
n0 = floor((1-over)*NN);
nsect=floor((N-NN)/(n0));
n1=1;

for dr=1:nsect
    idr=idr+1;
    decay=channel_3(ii,:);

%test=abs(decay+dr);[Y1,I1] = min(test);
%test=abs(decay+dr+10);[Y1,I2] = min(test)
I1=n1;
I2=(n1+NN-1);
sig=channel_3(ii,I1:I2)';
l_reg=(I2-I1+1);
m=30/l_reg;%start at estimated grad 
x=(I1:I2)';
A = [ones(l_reg,1),x];

%c = pinv(A)*sig;
%line=c(1)+c(2)*(1:N);
%line2=c(1)+c(2)*(1:length(impuls));
%line_s(ii,idr,:)=line;
%c_s(ii,idr,1:2)=c;
%I1_s(ii,idr)=I1;
%I2_s(ii,idr)=I2;

%figure
%plot(line)
n1=n1+n0;
end

%    hold on

%plot(t,10*log10(channel_2(ii,:)))
end


n1=1;
idr=0;
    %lin(1:N)=0;
if over ==0
for dr=1:nsect
    idr=idr+1;
    I1=n1;
    I2=(n1+NN-1);

 %   figure
  %  hist(c_s(:,idr,2))
   % I1=I1_s(1,idr);
   % I2=I2_s(1,idr);
    l_reg=(I2-I1+1);    
    decay=channel_3(:,I1:I2);
    [cc II]=min(sum(channel_3(:,I1:I2)'.^2));

      % min(c_s(:,idr,2));
  %  c1=(c_s(II,idr,1));
  %  cc=c_s(II,idr,2);
    %line2(idr,1:N)=nan;

    if idr==1
        %line2(idr,1:N)=c1+cc*(1:N);
        %lin(I1:I2)=line2(idr,I1:I2);
        %store=lin(I2);
          lin(I1:I2)=channel_3(II,I1:I2);
        
    end
    if idr>1
        %line2(idr,1:N)=c1+cc*(1:N);%
        %m=-1/(l_reg-1);
        %y=-m*(0:(l_reg-1))+1;
        %lin(I1:(I2))=line2(idr,I1:(I2));%+line2(idr,I1:(I2));
        
        %store+line2(idr,I1:(I2))-line2(idr,I1);
        %store=lin(I2); 
        lin(I1:I2)=channel_3(II,I1:I2);
    end
    
            
        
    n1=n1+n0;    
end
end


idr=0;

if over>0
for dr=1:nsect
    idr=idr+1;
 %   figure
  %  hist(c_s(:,idr,2))
      I1=n1;
    I2=(n1+NN-1);

    %I1=I1_s(1,idr);
  %  I2=I2_s(1,idr);
    l_reg=(I2-I1+1);    
    decay=channel_3(:,I1:I2);
    [cc II]=min(sum(channel_3(:,I1:I2)'.^2));
      % min(c_s(:,idr,2));
%    c1=(c_s(II,idr,1));
   % cc=c_s(II,idr,2);

    %line2(idr,1:N)=nan;

    if idr==1
        %line2(idr,1:N)=c1+cc*(1:N);
        %lin(I1:I2)=line2(idr,I1:I2);
        %store=lin(I2);
        lin(I1:I2)=channel_3(II,I1:I2);
        last=squeeze(channel_3(II,1:end));
    end
    if idr>1
        %line2(idr,1:N)=c1+cc*(1:N);%
        %m=-1/(l_reg-1);
        %y=-m*(0:(l_reg-1))+1;
        %lin(I1:(I2))=line2(idr,I1:(I2));%+line2(idr,I1:(I2));
        
        %store+line2(idr,I1:(I2))-line2(idr,I1);
        %store=lin(I2); 
        window=hanning(l_reg);
        win=window;
        win(l_reg/2:l_reg)=1;
        win2=abs((win)-1);
        %figure
        %plot(win2,'r')
        %hold on
        %plot(win)
        %pause(1)
           current=channel_3(II,I1:I2);
        lin(I1:I2)=last(I1:I2).*win2'+current.*win';
        last=channel_3(II,:);;
    end
    
            
        
       n1=n1+n0; 
end
end
%channel_final=min(channel_2,[],1);
%channel_final_back(length(channel_final):-1:1)= cumsum(channel_final((length(channel_final(:)):-1:1)).^2);
clear channel_final_back
%channel_final_back=lin;

line_lin=lin;%10.^(lin/10);
%plot(10*log10(channel_final_back/max(channel_final_back)),'red')
if (plot_yn==1)
    channel_final_back(length(line_lin):-1:1)= cumsum(line_lin((length(line_lin(:)):-1:1)).^2);
channel_final_back=channel_final_back/max(abs(channel_final_back));

plot(10*log10(channel_final_back),'red')
plot(10*log10(channel_fil1),'g')
ylim([-45 0])
end
%t1=(0:(length(channel_fil1)-1))/Fs
%plot(t1,10*log10(channel_fil1/max(channel_fil1)),'green')
%ylim([-60 0])
%xlim([0 2])

%figure

%hold on

%plot(lin','LineWidth',1)
%hold on%
%plot(10*log10(channel_final_back))
%channel_final=((channel_final));
channel_final=line_lin;
%xlim([0 length(channel_fil)])
%figure
%plot(y1)
%hold on
%plot(y2,'r')

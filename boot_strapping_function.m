function [T25_boot EDT_boot C80_boot centre_boot Y_logs Y_log]=boot_strapping_function(channel_opt_store,Fs)

%  xxx=-100:0.5:0;
%  xxx_fine=-100:0.01:0;
xxx=0.00000000000000000000001:0.001:1;
xxx_fine=0.00000000000000000000001:0.0001:1;
xxx_fine=0.00000000000000000000001:0.0005:1;
%logspace(1,10^(-100/20));
position=0;    
channel_final_store=channel_opt_store;
ll=size(channel_final_store)
for iii=1:ll(1)
%     str=sprintf('RNCM/claps%g.wav',iii)
%     if (iii==1)    str=sprintf('bradford/snare_1a.wav',iii);end
%     if (iii==2)    str=sprintf('bradford/snare_1b.wav',iii);end
%     if (iii==3)    str=sprintf('bradford/snare_1c.wav',iii);end
%     if (iii==4)    str=sprintf('bradford/snare_1d.wav',iii);end
%     str
%     [y Fs2]=wavread(str);
%     
% 
%     [m I]=max(abs(y(:,1)));y=y(I:end,1);
% 
%     y = resample(y,Fs,Fs2)';
%      Nf=6;
%      Fc=1000
%     wc=[Fc/sqrt(2) (Fc/sqrt(2))*2]/(Fs/2);
%     [B1,A1] = butter(Nf,wc);
%     clear channel_fil2
%     y=filter(B1,A1,y);
% if (iii==1 | iii==3)
%    
%     y_win=abs(y(1:1.4*Fs));
% end
% if (iii==2 | iii==4)
%    
%     y_win=abs(y(1:0.6*Fs));
% end
%                 [a,b,alpha]=MLE_3_function(y_win,Fs);
% %            toc
%              N=3*Fs;
%             channel_est=(alpha*a.^(0:(N-1))+(1-alpha)*b.^(0:(N-1)));
%             channel_final_store(iii,:)=channel_est;
            channel_next=channel_final_store(iii,:);
             channel_next1(length(channel_next):-1:1) = cumsum(channel_next(length(channel_next):-1:1).^2);
            Y_log(iii,:)=10*log10(channel_next1/max(abs(channel_next1)));
            plot(Y_log(iii,:))
            pause(0.1)
end
%%
position=1
lim(position,1:2)=[1;1];
gaps=10;
len_est=3;
% dist_decay_Y=zeros(100,length(xxx_fine));
for tt=gaps:gaps:len_est*Fs
 tt/(len_est*Fs)
    position=position+1;
   time(position)=(tt-1)/Fs;


    III=find((Y_log(:,tt)<0));
    %channel_final_store(III,tt)=Y_log(III,tt);
   % dist_decay(position,:)=hist(Y_log(III,tt),xxx);    
%     [parmhat, parmci] = evfit(dist_decay(position,:));
     %Y(position,:) = evpdf(xxx, parmhat(1), parmhat(2));
% save brad_dist_decay_
    %[parmhat, parmci] = lognfit(-Y_log(III,tt));
Yfit = ksdensity(channel_final_store(III,tt),xxx);
Ys(position,:)=Yfit;%/sum(Yfit);
%    cdf=cumsum(Y(position,:)/sum(Y(position,:)));
%    [m I05]=min(abs(0.05-cdf));
%    [m I95]=min(abs(0.95-cdf));
%    lim(position,1:2)=[xxx(I05);xxx(I95) ];
channel_final_store(III,:)';
[bootstat,bootsam] = bootstrp(500, 'median', channel_final_store(III,tt));  
%[bootstat,bootsam] = bootstrp(500, 'median', Y_log(III,tt));  
 Y = ksdensity(bootstat,xxx_fine);
    dist_decay_Y(position,:)=Y;    

%Y=hist(bootstat,xxx_fine);
%  plot(Y)
%  pause(0.01)
cdfm=cumsum(Y/sum(Y));
plot(cdfm)
[m I05]=min(abs(0.025-cdfm));
[m I95]=min(abs(0.975-cdfm));
lim(position,1:2)=[xxx_fine(I05);xxx_fine(I95) ];lim(position,1:2);

%     Y(position,:) = lognpdf(-xxx, parmhat(1), parmhat(2));
%     bar(xxx,dist_decay(position,:))
% 
% %     plot(xxx,Y(position,:),'r')
% 
%     pause(0.5)
%[h,p,ci] = ttest(Y_log(III,tt)-median(Y_log(III,tt)),0)
end
save rncm_dist_decay dist_decay_Y xxx_fine
%%
% [bootstat,bootsam] = bootstrp(100, 'median', Y_log(III,tt));
% m=hist(bootstat,xxx);
% 
% cdfm=cumsum(m/sum(m));
% [m I05]=min(abs(0.05-cdfm));
% [m I95]=min(abs(0.95-cdfm));
% lim(position,1:2)=[xxx(I05);xxx(I95) ];lim(position,1:2)

channel_ests=median(channel_final_store(:,:));
 channel_est(length(channel_ests):-1:1) = cumsum(channel_ests(length(channel_ests):-1:1).^2);
Y_log=10*log10(channel_est/max(abs(channel_est)));
clear Y_logs
 for pp=1:2
     channel_est=lim(:,pp);
     channel_est(length(channel_est):-1:1) = cumsum(channel_est(length(channel_est):-1:1).^2);
    Y_logs(:,pp)=10*log10(channel_est/max(abs(channel_est)));

 end
Ts_bootstraps=gaps*(1/Fs);Fs_bootstraps=1/Ts_bootstraps;

[T25_boot(2) EDT_boot(2) C80_boot(2) C50_boot(2) centre_boot(2) D_boot(2) ]=Room_acoustic_params_centre_ldr(channel_ests,Fs,25);

[T25_boot(1) EDT_boot(1) C80_boot(1) C50_boot(1) centre_boot(1) D_boot(1) ]=Room_acoustic_params_centre_ldr(lim(:,1)',Fs_bootstraps,25);
  [T25_boot(3) EDT_boot(3) C80_boot(3) C50_boot(3) centre_boot(3) D_boot(3) ]=Room_acoustic_params_centre_ldr(lim(:,2)',Fs_bootstraps,25);      
T25_boot'
EDT_boot'
C80_boot'
centre_boot'
D_boot'
figure
plot(time,Y_logs)
hold on
plot((0:(length(Y_log)-1))/Fs,Y_log,'r')
   ylim([-35 0])
     xlim([0 2])
     xlabel('Time (s)')
     ylabel('Level(dB)')
     legend('Lower 95% confidence limit','Upper 95% confidence limit','Median decay curve')
%%
% tt=500
%     III=find((Y_log(:,tt)<0));
% data=Y_log(III,tt);

% 
% [XI,YI] = meshgrid(time,10*log10(xxx.^2));
% ZI = griddata(time,10*log10(xxx.^2),Ys',XI,YI);
% figure
% 
% %surface(Y')
% 
% 
% % plot(Y(position,:))
% surf(XI,YI,ZI,'EdgeColor','none')
% view(2)
% shading('interp')
% %caxis([-0.001 1 ])
% zlim([-1 100])
% ylim([-60 0])
   clear all
close all
% matlabpool local 8
%  load Robs_files_curtains_close
Fs2=3000;
Fcm=[63    125 250   500 1000 2000 4000  8000];
Fs3=[3000 3000 3000 3000 3000 6000 12000 24000];
% Fcm=[63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 ];
% for i=1:length(Fcm)
%     if Fcm(i)<1001 Fs3(i)=3000; end
%     if Fcm(i)>1001 & Fcm(i)<2001 Fs3(i)=6000; end
%      if Fcm(i)>2001 & Fcm(i)<4001 Fs3(i)=12000; end
%     if Fcm(i)>4001 Fs3(i)=24000; end
% end
% imp=1
% % Fs3=[3000 3000 3000 3000 3000 6000 12000 24000];
% %%
% for maind=1:3
% if maind==1 direcmain=sprintf('DS_English_E8/');end
% % if maind==2 direcmain=sprintf('../Addenbrookes/Ward N3/');end
% if maind==2 direcmain=sprintf('DS_Maths_M7/');end
% if maind==3 direcmain=sprintf('DS_Science_S6/');end
% 
% files22=dir(direcmain);
% for lesson=1%:length(files22)
%     direc1=sprintf('%s',direcmain)
% % direc=sprintf('../Addenbrookes/Ward D8/');
% files11=dir(direc1);
% % for day=3:length(files11)
% 
% 
% 
% impulse_esti=zeros(8,12000);
% % direc=sprintf('%s/%s/',direc1,files11(day).name);
% files=dir(direc1);
% LL=length(files)
% ii=0;
% clear time_file
% for i=3:LL
%     ii=ii+1;
%     n=files(i).date;
%         sprintf('%s+%s*60',n(19:20),n(16:17));
% 
%     time_file(ii)=eval(sprintf('%s+%s*60+%s*60*60;',n(19:20),n(16:17),n(13:14)));
% end
% [Y,Ifiles] = sort(time_file);
%     %filename=sprintf('English_091110_BDB_11025.wav'); str=sprintf('speech')
%     %filename=sprintf('Maths_091109_BDB_11025.wav'); str=sprintf('speech')
%      % filename=sprintf('Science_091111_BDB_11025.wav'); str=sprintf('speech')
% 
% %     filename=sprintf('English_091110_BDB_24000.wav'); str=sprintf('speech')
%     %filename=sprintf('Maths_091109_BDB_24000.wav'); n(1)str=sprintf('speech')
%       %filename=sprintf('Science_091111_BDB_24000.wav'); str=sprintf('speech')
% %%    
% 
% for Fci=2:7;%3:(length(Fcm)-1)
% a_store1=[]
%    b_store1=[]
%    alpha_store1=[]
%    len_store1=[]
%    DR1=[]
% for filenum=1:length(time_file)
%     [Fcm(Fci) 100*filenum/length(time_file)]
%   
% 
%     %filename=sprintf('English_091110_BDB_11025.wav'); str=sprintf('speech')
% %     filename=sprintf('Maths_091109_BDB_11025.wav'); str=sprintf('speech')
% clear filename
% 
% direc_store=direc1;
%         filename=sprintf('%s%s',direcmain,files(Ifiles(filenum)+2).name);
% 
%     %filename=sprintf('Science_091111_BDB_11025.wav'); str=sprintf('speech')
% % if Fci>5
% %     %filename=sprintf('English_091110_BDB_24000.wav'); str=sprintf('speech')
% %     filename=sprintf('Science_091120_Weydon2_12000.wav'); str=sprintf('speech')
% % %     filename=sprintf('Maths_091109_BDB_24000.wav'); str=sprintf('speech')
% %     %filename=sprintf('Science_091111_BDB_24000.wav'); str=sprintf('speech')
% % end
%   SIZ=wavread(filename,'size'); 
%     Fs2=Fs3(Fci);
% 
% 
%     split_=1;
%     toN=0
%     window_split=round(SIZ(1)/split_);
%     for split_wav=1:split_
%         split_wav/split_;
%         fromN=toN+1;
%         toN=toN+window_split;
%         if toN> SIZ  toN=SIZ(1);end
%     [music,Fsw,NBITSw]=wavread(filename,[fromN toN]);%input signal
%     
%     
%     music=music(:,1);
%     music = resample(music,Fs2,Fsw)';%resanmple to Fs2
% %     music=music(1:2*61*Fs2);
%     TT(filenum)=length(music)/Fs2/60;
% %     sec=round(TT/1)    
%     sec=1;
% 
%  len_sec=length(music);
%     Fc=Fcm(Fci);
%     Nf=5;
% %     F1=2^(-1/(2))*1000
% %     F2=2^(-1/(2))*1000
%      wc=[Fc/sqrt(2) (Fc/sqrt(2))*2]/(Fs2/2);
% %         wc=[Fc/(10^0.05) Fc*(10^0.05)]/(Fs2/2);
% %     1000*10^0.05
%     [B1,A1] = butter(Nf,wc);
%     music_fil=filter(B1,A1,music);clear music
% %%  Segmentation 
%      %%input is music_fil(NxM)  N is the number of sections / recordings M is
%      %%the length of the sections (4 mins works well) but dont have to be
%      %%equal     
% 
% %
% %%
% 
% 
% %%
% % music_fil1=music_fil(1:60*Fs2);
% 
% 
% %   [channel_final_25(filenum,:) channel_final_20(filenum,:)  impulse_est(filenum,:) start_end2(filenum,2,:) a_store(filenum,:) b_store(filenum,:) alpha_store(filenum,:) len_store(filenum,2,:) DR(filenum,:)]=MLE_decay_estimate(music_fil,sections,Fs2,Fc,imp);
% % % % % % % % % % % % % % % % 
% 
% [channel_final_25 channel_final_20  impulse_est start_end2 a_store b_store alpha_store len_store DR]=MLE_decay_estimatep(music_fil,1,Fs2,Fc,imp);
% % [channel_final_25 channel_final_20  impulse_est start_end2 a_store b_store alpha_store len_store DR]=MLE_decay_estimate(music_fil,sections,Fs2,Fc,imp);
% % [channel_final_25 channel_final_20  impulse_est start_end2 a_store b_store alpha_store len_store DR]=MLE_decay_estimate(music_fil,sections,Fs2,Fc,imp);
%    a_store1=[a_store1; a_store];
%    b_store1=[b_store1; b_store];
%    alpha_store1=[alpha_store1; alpha_store];
%    len_store1=[len_store1 ;len_store];
%    DR1=[DR1 DR];
%    
%    % 
% % size(channel_final_25)
% % size(len_store)
% 
%   % if length(impulse_est)==0
% %     impulse_est(1:4*Fs2)=nan;
% % end
% 
% %%%  Calculate parameter estiamtes from speech
% impulse_est=median(channel_final_25);
% [T25_est(imp) EDT_25_est(imp) C80_25_est(imp) C50_25_est(imp) centre_25_est(imp) D_25_est(imp) ]=Room_acoustic_params_centre_ldr(impulse_est,Fs2,25);
% [T25_est(imp) EDT_25_est(imp) C80_25_est(imp) C50_25_est(imp) centre_25_est(imp)];
% end
% %  N=2000;ib=2
% % impulse_est=(alpha_store(ib)*a_store(ib).^(0:(N-1))+(1-alpha_store(ib))*b_store(ib).^(0:(N-1)));
% % plot(10*log10((impulse_est.^2)))
% 
% % [T25_est(imp) EDT_25_est(imp) C80_25_est(imp) C50_25_est(imp) centre_25_est(imp) D_25_est(imp) ]=Room_acoustic_params_centre_ldr(impulse_est,Fs2,25);
% % [T25_est(imp) EDT_25_est(imp) C80_25_est(imp) C50_25_est(imp) centre_25_est(imp)]
% 
% % [T25(imp) EDT(imp) C80(imp) C50(imp) centre(imp)  ]
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% %  filename=sprintf('bookended_%g_arti_Fs2_%g_Fc_%g.mat',sections,Fs2,Fc)
% %  [T25_boot EDT_boot C80_boot centre_boot Y_logs Y_log]=boot_strapping_function(channel_final_25,Fs2);
% %  save(filename,'channel_final_20','channel_final_25','T25_est','EDT_25_est','C80_25_est','centre_25_est','impulse_est','start_end2','a_store', 'b_store' ,'alpha_store' ,'len_store','DR','Y_logs','Y_log')
% %%
% end
% 
% 
%    a_store1_cel{Fci}=a_store1;
%    b_store1_cel{Fci}=b_store1;
%    alpha_store1_cel{Fci}=alpha_store1;
%    len_store1_cel{Fci}=len_store1;
%    DR1_cel{Fci}=DR1; 
%    str=sprintf('b%sfc%g.mat',files(Ifiles(filenum)+2).name,Fcm(Fci))
%    save(str,'a_store1_cel','b_store1_cel','alpha_store1_cel','len_store1_cel','DR1_cel');
% 
% end
% 
% %  save Bedford_Bay1_06_00_10_00.matlesson
% 
% end 
% end
   clear all
   close all
   
   
   %%load('\\D6hrk04j\shareddocs\PK\schools_BMLERT\bDS_English_E8_L5resamp.wavfc1000.mat')

% break
%%  
%    ho=hi   
 for Fci=3:7
    
%    plot(-80:20:0,ho)
%    hold on
%    plot(-80:20:0,hc,'r')
% % 
%%
Fs2=3000;
Fcm=[63    125 250   500 1000 2000 4000  8000];
Fs3=[3000 3000 3000 3000 3000 6000 12000 24000];
% Fci=5
Fs2=Fs3(Fci);
% bDS_English_E8_L5resamp.wavfc250
% bDS_Science_S6_L5resamp.wavfc1000.wav
% bDS_Maths_M7_L5resamp.wavfc1000

 str=sprintf('ordered_bDS_Science_S6_L5resamp.wavfc%g.mat',Fcm(Fci))
         load(str,'a_store1_cel','b_store1_cel','alpha_store1_cel','len_store1_cel','DR1_cel','pos_store1_cell');
% str=sprintf('posDS_Maths_M7_L5resamp.wavfc%g.mat',Fcm(Fci))
%          load(str);

   DR1=DR1_cel{Fci};
   a_store1=a_store1_cel{Fci};
   b_store1=b_store1_cel{Fci};
   alpha_store1=alpha_store1_cel{Fci};
   len_store1=len_store1_cel{Fci};
   pos_store1=pos_store1_cell{Fci};
%    Rand_order=randperm(length(DR1));
%    
%      DR1=DR1(Rand_order);
%    a_store1=a_store1(Rand_order);
%    b_store1=b_store1(Rand_order);
%    alpha_store1=alpha_store1(Rand_order);
%    len_store1=len_store1(Rand_order);
   
   III=find(DR1<-25);
   

% st(DR1,-80:20:0)
if length(III)>12
%   if  direc(38:42)=='0600'
%       sLen
% sLen=;floor((sum(TT)*60)/(60*2));
%   pos_less_start=[ 1 4659 9831 13810 16760 inf]; Science!
for i=2:(length(pos_store1)-1)
    id(i)=0;
    if (pos_store1(i)<pos_store1(i-1))
id(i-1)=1;
    end
end
posI=find(id==1);
pos_less_start=[0 0 0 0  inf]'
pos_less_start(1)=1;pos_less_start(2:4)=(posI);
sLen=40
section=floor(length(III)/sLen)
clear channel_store_sec 
channel_store_sec=zeros(sLen,4*Fs3(Fci));
Y_s=channel_store_sec;

for lesson=1:4
close all
    III=find(DR1<-25); 
 i_s=   find(III>pos_less_start(lesson) & III<pos_less_start(lesson+1));

      III1=find(DR1(III(i_s))<-25); 
sLen=8
section=round(length(III1)/sLen)
for sec=1:sLen
  
     pos1=(sec-1)*section+1;
    pos2=(sec-1)*section+1+section-1;
    if sLen==sec 
        pos2=length(III1);
    end
  i=(pos1:pos2);  [pos1 pos2]
%  
%  for ii=1:5
%  if(III(pos1)>pos_less_start(ii))% & III(pos2)<pos_less_start(ii+1))
%     time_end(sec)=ii;
%  end
%  end

%  errorbar(1:5,medRT25_room(Fci,:),sdRT25_room(Fci,:),'x')
 
 
    channel_fil=zeros(3*Fs2,1);
    [channel_sec]=optimum_model_function_2(a_store1(III(i_s(i))),b_store1(III(i_s(i))),alpha_store1(III(i_s(i))),Fs2,DR1(III(i_s(i))),channel_fil,0);

%     Find -30dB point
channel_sec_log=10*log10(channel_sec.^2);
% [mm Iend]=min(abs(channel_sec_log-(-25)));
% channel_sec=channel_sec(1:Iend);

%     Find -30dB point
channel_sec_log=10*log10(channel_sec.^2);
[mm Iend]=min(abs(channel_sec_log-(-30)));
channel_sec=channel_sec(1:Iend);
channel_sec2=channel_sec;
Y_log=channel_sec_log;
%find point where sound has decayed 5dB
test=abs(Y_log+5);
[Y1,I1] = min(test);
%find point where sound has decayed 35dB
test=abs(Y_log+35);
[Y2,I2] = min(test);
%%%%%%%%%%least squared error T30
sig=Y_log(I1:I2)';
l_reg=(I2-I1+1);
% m=30/l_reg;%start at estimated grad 
x=(I1:I2)';
A = [ones(l_reg,1),x];
c = pinv(A)*sig;
line=(channel_sec_log(I2))+c(2)*((1:length(channel_sec))-I2);
line=10.^(line(:)/20);
% plot(20*log10(line))

plot(channel_sec_log,'r')
hold on
winL=500;
channel_sec2(I2:end)=line(I2:end)';
plot(20*log10(channel_sec2),'k')
 channel_sec=channel_sec2;



    [T25_sec(lesson,Fci,sec) EDT_sec(Fci,sec) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_sec,Fs2,25);
channel_store_sec(sec,1:length(channel_sec))=channel_sec;

% [T25_sec(sec) EDT_sec(sec) C80 C50 centre D ]

%              N=6*Fs2;
%             channel_est=(alpha_store1(III(i))*a_store1(III(i)).^(0:(N-1))+(1-alpha_store1(III(i)))*b_store1(III(i)).^(0:(N-1)));
%             clear Y

             Y(length(channel_sec):-1:1) = cumsum(channel_sec(length(channel_sec):-1:1).^2);
             Y=Y/max(abs(Y));
             Y_s(sec,1:length(Y))=Y;
             
%           
% %             DR(ii)=Y_log(length(y_re));
%             [T25_e(i) EDT_e(i) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_est,Fs2,25);
% 
    

end
%%
t=((1:length(Y_s))-1)/Fs2;
figure
plot(t,10*log10(Y_s'))
ylim([-35 0])
channel_final=median(channel_store_sec);
% impulse_esti(Fci,:)=median(channel_store_sec);
impulse_est_cell{Fci}=median(channel_store_sec);
%  for i=1:5
% %     impulse_est_cell{i}=impulse_esti(i,:);
% Y_logs_cell{i}=Y_logs(i,:,1:2);
% Y_log_cell{i}=Y_log(i,:);
%  end

%  [T25_boot(Fci,:) EDT_boot(Fci,:) C80_boot centre_boot Y_logs_cell{Fci} Y_log_cell{Fci}]=boot_strapping_function(channel_store_sec,Fs2);
% [channel_final]=optimum_model_function_2(a_store1(III(:)),b_store1(III(:)),alpha_store1(III(:)),Fs2,DR1(III(:)),channel_fil,1)
[T25_e(Fci) EDT_e(Fci) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_final,Fs2,25);
% close all
% [T25_e EDT_e C80 C50 centre D ]=Room_acoustic_params_centre_ldr_plot(channel_final,Fs2,25);
% uiopen('C:\Documents and Settings\AGS056\My Documents\Thesis\schools_recordings\impulse_english.wav',1)
%  %%
% for Fci=5
%     
% wc=[Fcm(Fci)/sqrt(2) (Fcm(Fci)/sqrt(2))*2]/(Fs2/2);
% %         wc=[Fc/(10^0.05) Fc*(10^0.05)]/(Fs2/2);
% %     1000*10^0.05
% Nf=5
%     [B1,A1] = butter(Nf,wc);
%     dataf=filtfilt(B1,A1,data);
% [mm datafI]=max(abs(dataf));
% dataf=dataf(datafI:end);
% [T25_i(Fci) EDT_i(Fci) C80i C50i centrei Di ]=Room_acoustic_params_centre_ldr_plot(dataf',fs,15);
% end


end
end
%  for roomi=1:5
%      IIi=find(time_end==roomi);T25_sec(Fci,IIi)'
%      RT25_room{Fci,roomi}=T25_sec(Fci,IIi);
% 
%  end


 end

% end
% % 06.00-10.00
% %  save Bedford_Bay1_06_00_10_00.mat
% % str=sprintf('%s/day_%g.mat',direc_store,day);
% %  save(str)
% end
% end
%  end
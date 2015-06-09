% Blind estimation of reverberation time, paul kendrick, University of
% Salford, Acoustics research centre
clear all
close all
%Octave bands
Fcm=[63  125 250   500 1000 2000 4000  8000];
%Resample to each Fs for each octave band
Fs3=[3000 3000 3000 3000 3000 6000 12000 24000];

%%  this section loops over a wav file and extracts all decay phases,
%  fits a Maximum likelihood model to each decay phase, computes the
%  dynmaic range, and collects all the results for the following wav file
filename=sprintf('originals/DS_English_E8_L1.WAV');
SIZ=wavread(filename,'size');

for Fci=5%:8
    % clear Maximum likelihood parameter
    a_store1=[];
    b_store1=[];
    alpha_store1=[];
    len_store1=[];
    DR1=[];
    start_end1=[];
    
    Fs2=Fs3(Fci);  %sampling freqeuency to use for this octave band
    
    % reads sections of wav file at a time to save memory, split_ is the
    % number of sections
    split_=20;    toN=0;
    window_split=round(SIZ(1)/split_);
    for split_wav=1:split_
        %% read bit of wav file in from fromN to toN, then filter
        fromN=toN+1;
        toN=toN+window_split;
        if toN> SIZ  toN=SIZ(1);end
        [music,Fsw,NBITSw]=wavread(filename,[fromN toN]);%input signal
        music=music(:,1);
        music = resample(music,Fs2,Fsw)';%resanmple to Fs2
        Fc=Fcm(Fci);
        Nf=5;
        wc=[Fc/sqrt(2) (Fc/sqrt(2))*2]/(Fs2/2);
        [B1,A1] = butter(Nf,wc);
        music_fil=filter(B1,A1,music);clear music
        %% This function runs the polyfit algorithm, finds decay sphases 
        %and extracts the ML parameters for the decay model
        [start_end_store a_store b_store alpha_store len_store DR ]=MLE_decay_estimatep(music_fil,1,Fs2,Fc,1);
 
        %concatonate all results
        a_store1=[a_store1 a_store];
        b_store1=[b_store1 b_store];
        alpha_store1=[alpha_store1 alpha_store];
        len_store1=[len_store1 len_store];
        DR1=[DR1 DR];
        start_end1=[start_end1 start_end_store+fromN-1];
    end
    %  filename=sprintf('bookended_%g_arti_Fs2_%g_Fc_%g.mat',sections,Fs2,Fc)
    %  [T25_boot EDT_boot C80_boot centre_boot Y_logs Y_log]=boot_strapping_function(channel_final_25,Fs2);
    %  save(filename,'channel_final_20','channel_final_25','T25_est','EDT_25_est','C80_25_est','centre_25_est','impulse_est','start_end2','a_store', 'b_store' ,'alpha_store' ,'len_store','DR','Y_logs','Y_log')
    a_store1_cel{Fci}=a_store1;
    b_store1_cel{Fci}=b_store1;
    alpha_store1_cel{Fci}=alpha_store1;
    len_store1_cel{Fci}=len_store1;
    start_end1_cel{Fci}=start_end1;
    DR1_cel{Fci}=DR1;
    
    % str=sprintf('b%sfc%g.mat',files(Ifiles(filenum)+2).name,Fcm(Fci))
    % save(str,'a_store1_cel','b_store1_cel','alpha_store1_cel','len_store1_cel','DR1_cel');
    
end

%  save Bedford_Bay1_06_00_10_00.matlesson
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
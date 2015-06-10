% Blind estimation of reverberation time, Paul kendrick, University of
% Salford, Acoustics research centre, p.kendrick@salford.ac.uk
% This program can compute the reverberation time from passivley recieved
% signals, such as speech and music. The operation is as follows:
%   Description of algorithm
%   1) Choose wav file, this needs to be a recording of speech or music,
%   for speech at least an hour, for music perhaps longer
%   2)  read in section of wav file, and filter into octave band
%   3)     Maximum likehood estimation of decays; 
%          a) Polyfit algorithm, this uses 0.5 s windows on the envelope of
%           the signal, fitting a 1st order poly nomial to the log
%           envelope, the gradient of each fit is extracted, windows are
%           moved along, using a very high overlap (0.002s), and regions
%           which are continuously decaying identified from the gradient
%           (gradient equivlent for RT==100s used as maximum gradient) -
%           these decay phases are then fine tuned (start maximum, end min
%           of LP filtered envelope)
%          b) Maximum Likelihood fit to every decay phase, a model of sound
%          decay in a room is where the envelope, alpha*a.^I+(1-alpha)*b.^I,
%          modulates gaussian white noise.  This is a
%          model of non-diffuse sound decay, utilising two exponential
%          decays, with decay rates determined by a and b, added together
%          using a convex sum where the balence between the two is
%          controlled by the parameter alpha, all parmaters are constrained
%          between 0 & 1. The fit proceedure first computes the likelihood
%          over a coarse grid of values for a and b (alpha optimsed for
%          each grid point), then performs a fine search using the corase
%          result as a starting point, the function fmincon is used.  This
%          is contrained minimisation.
%          c)   Once all ML parameters are computed, the dynamic range of
%          each decayphase is computed from the ML decay curve.
%   4)     Postprocess results, for every octave band, a framework for
%   estimationg the RT is as follows, find the length of signal required to
%   ensure that at least 40 decay phases with at least 25 dB of dynamic
%   range are present, (this is just a rule of thumb), I found for speech
%   3-5 (segLen is this parameter - in secs) mins works.  In each 3-5 min 
%  section, compute the RT from  the 
%  decay phases with at least 25dB dynamic range, select the minimum RT as
%  the estimate for that section.  Compute the minimum RT for multiple
%  sections.  Average over all sections provides the blind RT estimate, the
%  95% CL from the standard error indicates the confidence in that result.
%   

clear all
close all
%Octave bands
Fcm=[63  125 250   500 1000 2000 4000  8000];
%Resample to each Fs for each octave band
Fs3=[3000 3000 3000 3000 3000 6000 12000 24000];

%%  this section loops over a wav file and extracts all decay phases,
%  fits a Maximum likelihood model to each decay phase, computes the
%  dynmaic range, and collects all the results for the following wav file
filename=sprintf('originals/DS_Science_S6_L1.WAV');
SIZ=wavread(filename,'size');

for Fci=1:8
    
    % clear Maximum likelihood parameters
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
    splitFs2=1;
    for split_wav=1:split_
        split_wav
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
        start_end1=[start_end1 start_end_store+splitFs2-1];
        splitFs2=length(music_fil)+splitFs2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute RT for room
% choose Length of segments
segLen=3*60;  %% three mins (in secs)

clear fittedresults
for Fci=1:5
    
    Fs2=Fs3(Fci);
    segLenN=segLen*Fs2;
    
    DR1=DR1_cel{Fci};
    a_store1=a_store1_cel{Fci};
    b_store1=b_store1_cel{Fci};
    alpha_store1=alpha_store1_cel{Fci};
    len_store1=len_store1_cel{Fci};
    pos_store1=start_end1_cel{Fci};% pos_store1_cell{Fci};
    start_end1=pos_store1;
    lenFile=max(max(start_end1));
    lenFileS=max(max(start_end1))/Fs2;
    
    Nseg=ceil(lenFile/segLenN);
    
    % find all decay phases with DR better than 25dB
    III_=find(DR1<-25);
    %%      
    clear T25_sec NdecayPhases Y_s
    channel_store_sec=zeros(Nseg,4*Fs3(Fci));
    channel_fil=zeros(5*Fs2,1);
    for segi=1:Nseg
        
        posN_from=(segi-1)*segLenN+1;
        posN_to=(segi)*segLenN;
        if segi==Nseg
            posN_to=lenFile;
        end
        
        III=III_(start_end1(1,III_)>=posN_from &  start_end1(2,III_)<posN_to);
        
        if length(III>1)
        [channel_sec]=optimum_model_function_2(a_store1(III),b_store1(III),alpha_store1(III),Fs2,DR1(III),channel_fil,0);
        channel_sec3=channel_sec;
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
        %       plot(20*log10(line))
        %       plot(channel_sec_log,'r')
        %       hold on
        %       winL=500;
        %       channel_sec2(I2:end)=line(I2:end)';
        %       plot(20*log10(channel_sec2),'k')
        %       channel_sec=channel_sec2;

        [T25_sec(segi) EDT_sec(Fci,segi) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_sec3,Fs2,25);
        channel_store_sec(segi,1:length(channel_sec3))=channel_sec3;
        else
            channel_store_sec(segi,1:length(channel_sec3))=channel_sec3.*nan;
            T25_sec(segi)=nan;
        end
        
        %         Y(length(channel_sec3):-1:1) = cumsum(channel_sec3(length(channel_sec3):-1:1).^2);
        %         Y=Y/max(abs(Y));
        %         Y_s(segi,1:length(Y))=Y;
        NdecayPhases(segi)=length(III);
    end
    
    % Collect together results
    
    fittedresults{Fci}.DR1=DR1;
    fittedresults{Fci}.a_store1=a_store1;
    fittedresults{Fci}.b_store1=b_store1;
    fittedresults{Fci}.alpha_store1=alpha_store1;
    fittedresults{Fci}.len_store1=len_store1;
    fittedresults{Fci}.pos_store1=pos_store1;% pos_store1_cell{Fci};
    fittedresults{Fci}.T25_sec=T25_sec;
    fittedresults{Fci}.NdecayPhases=NdecayPhases;
    
    %     t=((1:length(Y_s))-1)/Fs2;
    %     figure
    %     plot(t,10*log10(Y_s'))
    %     ylim([-35 0])
    %     channel_final=median(channel_store_sec);
    %     hold on
    %     t=((1:length(channel_final))-1)/Fs2;
    %     plot(t,20*log10(channel_final'),'linewidth',2)
    
    fittedresults{Fci}.channel_meadian=median(channel_store_sec);
    
    %     [T25_boot(Fci,:) EDT_boot(Fci,:) C80_boot centre_boot Y_logs_cell{Fci} Y_log_cell{Fci}]=boot_strapping_function(channel_store_sec,Fs2);
    [T25_e_decaycurve(Fci) EDT_e(Fci) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_final,Fs2,25);

    fittedresults{Fci}.T25_est = nanmean(T25_sec);
    fittedresults{Fci}.T25_cl = nanstd(T25_sec)/sqrt(sum(isnan(T25_sec)==0))*1.96;
    fittedresults{Fci}.T25_N=sum(isnan(T25_sec)==0);
   
    [fittedresults{Fci}.T25_est fittedresults{Fci}.T25_cl]
   
end

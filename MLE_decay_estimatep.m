function [start_end_store a_store b_store alpha_store len_store DR ]=MLE_decay_estimate(music,Sections,Fs,Fc,name)
% channel_final_25 channel_final_20 channel_est channel_opt_store

%%%%%%%%%%Full Automatic frameswork for MLE
channel_final_25=[];
channel_final_20=[];
channel_est=[];
start_end_store=[];
a_store =[];
b_store =[];
alpha_store =[];
len_store =[];
DR =[];
channel_opt_store=[];
a_storet=[];
b_storet=[];
alpha_storet=[];
len_storet=[];

start_end2a=[];
start_end2b=[];
start_end_storea=[];
start_end_storeb=[];

%%%%%%%%Stage 1 segment signal
% Nf=6;
% wc=[Fc/sqrt(2) (Fc/sqrt(2))*2]/(Fs/2);
% [B1,A1] = butter(Nf,wc);

input=music;
no_picks=100000;
y=input';
% compute Low pass filtered envelope
env=env_detect2(Fs,80,y);

% poly fit algorithm, finds regions which are contiuously decaying
wind=0.5;
maxRT=100;
step_size_s=0.00200;
start_end_total=Choose_signal(log10(env),Fs,no_picks,wind,maxRT,step_size_s);
no_picks_total=length(start_end_total);
t=(0:(length(env)-1))/Fs;
start_end=start_end_total;
si_no_picks=size(start_end);
no_picks=si_no_picks(2);
env2=env;%env_detect2_windowed(y,Fs,80,y,40);
env=abs(y);

%%%%%%%%%%%%%Stage 2 perform ML of segements
% te=sprintf('Please wait...performing on MLE Section %g of %g',section,Sections)
% h = waitbar(0,te);
%=2;
start_end2=zeros(2,no_picks);
start_end2a=zeros(1,no_picks);
start_end2b=zeros(1,no_picks);
start_end_storea=zeros(1,no_picks);
start_end_storeb=zeros(1,no_picks);
DR=zeros(1,no_picks);
T25=zeros(1,no_picks);
EDT=zeros(1,no_picks);
%  while ii<no_picks%
% waitbar(ii/no_picks)
N=6*Fs;
Y=zeros(1,N);

%% ii loops over the decay phases found by the polyfit algorithm
parfor ii=1:no_picks
    
    % This extracts the selected decay phase
    N1=abs(start_end(1,ii)-start_end(2,ii));
    y_re1=y(start_end(1,ii):(start_end(2,ii)-1));
    env_small=env(start_end(1,ii):(start_end(2,ii)-1));
    env_small2=env2(start_end(1,ii):(start_end(2,ii)-1));

    %fine tune the decay phases, i.e. find the max, thats the start, and
    %find the min of the lp filtered envelope, thats the end
    [m,Iend]=min(env_small((N1-(wind/2)*Fs):N1));Iend=Iend+(N1-(wind/2)*Fs)-1;y_re=y_re1(1:Iend);Nmin=length(y_re);
    [m,Iwin]=max(env_small(1:round((wind/2)*Fs)));y_re=y_re((Iwin):end);N=length(y_re);
    start_end2a(ii)=start_end(1,ii)+Iwin-1;
    start_end2b(ii)=start_end(1,ii)+Iwin-1+N-1;
    y_re=y_re/max(abs(y_re));
    y_win=y_re;
    
    % store the fine tuned start and end locations
    start_end_storea(ii)=start_end(1,ii)+Iwin-1;
    start_end_storeb(ii)=start_end(1,ii)+Iwin-1+N-1;
    
    % maximum likelihood fit of decay model to data
    [a,b,alpha]=MLE_3_function(y_win,Fs);

    a_storet(ii)=a;
    b_storet(ii)=b;
    alpha_storet(ii)=alpha;
    len_storet(ii)=length(y_re);
    
    % compute the decay curve form the ML fitted model
    N1=6*Fs;
    channel_est=(alpha*a.^(0:(N1-1))+(1-alpha)*b.^(0:(N1-1)));
    Y=zeros(N1,1);
    Y(length(channel_est):-1:1) = cumsum(channel_est(length(channel_est):-1:1).^2);
    Y_log=10*log10(Y/max(abs(Y)));
%     plot(Y_log)
    %plot(y_re);hold on;plot(channel_est,'r')
    %Dynamic range of decay phase
    DR(ii)=Y_log(length(y_re));
    [T25(ii) EDT(ii) C80 C50 centre D ]=Room_acoustic_params_centre_ldr(channel_est,Fs,25);
% soundsc(y_re,Fs)
end
a_store=a_storet;
b_store=b_storet;
alpha_store=alpha_storet;
len_store=len_storet;

start_end_store(1,1:no_picks)=start_end_storea;
start_end_store(2,1:no_picks)=start_end_storeb;


%  hold off
% (length(len_store(:,section)))
%              close(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Stage 3 estimate Decay curve from sections
% clear T35_est_sec T25_est_sec T20_est_sec T15_est_sec T10_est_sec EDT_5_est_sec
% T35_est_sec=[];
% T25_est_sec=[];
% T20_est_sec=[];
% T15_est_sec=[];
% clear EDT_35 EDT_25 EDT_20 EDT_15 EDT_10
% clear corr_rev35 corr_rev25 corr_rev20 corr_rev15 corr_rev10  T30_psd_store
% test=length(find(len_store(:,section)>0));
%
%              secs=find(abs(DR(1:no_picks))'>9 & len_store(1:no_picks,section)>(0.0*Fs) );test=length(secs);
%             cc=9;
%             if test==0
%                 secs=find(abs(DR(1:no_picks))'>0 & len_store(1:no_picks,section)>(0.0*Fs) );test=length(secs);
%                 cc=4;
%
%             end%%%% if there are no sections then reduce the DR constrinat
%
%             for i=1:test;
%                 ig=secs(i);
%
%                 if (abs(DR(ig))>35); rtdr=35;
%                     elseif ((abs(DR(ig))<=35 & abs(DR(ig))>25)); rtdr=25;
%                     elseif(abs(DR(ig))<=25 & abs(DR(ig))>20); rtdr=20;
%                     elseif(abs(DR(ig))<=20 & abs(DR(ig))>15); rtdr=15;
%                     elseif(abs(DR(ig))<=15 & abs(DR(ig))>10); rtdr=10;
%                     elseif(abs(DR(ig))<=10 ); rtdr=cc;
%                 end
%                 N=6*Fs;
%                 channel_next=(alpha_store(ig,section)*a_store(ig,section).^(0:(N-1))+(1-alpha_store(ig,section))*b_store(ig,section).^(0:(N-1)));
%                 dr_est(i)=rtdr;
%                 dr_real(i)=DR(ig);
%                 T35_est_sec(i)=0;T25_est_sec(i)=0;T20_est_sec(i)=0;T15_est_sec(i)=0;T10_est_sec(i)=0;
%                 EDT_35(i)=0;EDT_25(i)=0;EDT_20(i)=0;EDT_15(i)=0;EDT_10(i)=0;
%                 corr_rev35(i)=0;corr_rev25(i)=0;corr_rev20(i)=0;corr_rev15(i)=0;corr_rev10(i)=0;
%               if(dr_est(i)==35) [T35_est_sec(i) EDT_35(i) C80_35(i) C50_35(i) centre_35(i) D_35(i)]=Room_acoustic_params_centre_ldr(channel_next,Fs,35);
%               end
%                if(dr_est(i)>=25)[T25_est_sec(i) EDT_25(i) C80_25(i) C50_25(i) centre_25(i) D_25(i)]=Room_acoustic_params_centre_ldr(channel_next,Fs,25);
%
%                end
%                if(dr_est(i)>=20)[T20_est_sec(i) EDT_20(i) C80_20(i) C50_20(i) centre_20(i) D_20(i) ]=Room_acoustic_params_centre_ldr(channel_next,Fs,20);
%
%                end
%                if(dr_est(i)>=15)[T15_est_sec(i) EDT_15(i) C80_15(i) C50_15(i) centre_15(i) D_15(i) ]=Room_acoustic_params_centre_ldr(channel_next,Fs,15);
%
%                end
%                if(dr_est(i)>=10)[T10_est_sec(i) EDT_10(i) C80_10(i) C50_10(i) centre_10(i) D_10(i) ]=Room_acoustic_params_centre_ldr(channel_next,Fs,10);
%
%                end
%
%                 T35_est_store(i)=T35_est_sec(i);EDT35_est_store(i)=EDT_35(i);
%                 T25_est_store(i)=T25_est_sec(i);EDT25_est_store(i)=EDT_25(i);
%                 T20_est_store(i)=T20_est_sec(i);EDT20_est_store(i)=EDT_20(i);
%                 T15_est_store(i)=T15_est_sec(i);EDT15_est_store(i)=EDT_15(i);
%                 T10_est_store(i)=T10_est_sec(i);EDT10_est_store(i)=EDT_10(i);
%
%
%         end
%
%             channel_fil=1;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              I_dr=find(T35_est_sec>0);if size(T35_est_sec(I_dr)>0) [T35_est(section) I1]=min(T35_est_sec(I_dr));%start_end_35(1:2)=start_end(:,secs(I_dr(I1)));%%T35(section)=T35rt(I1);
%               C80_35_est(section)=min(C80_35(I_dr)); C80_35_estm(section)=median(C80_35(I_dr));C80_35_estmax(section)=max(C80_35(I_dr));
%               centre_35_est(section)=min(centre_35(I_dr));centre_35_estm(section)=median(centre_35(I_dr));
%               [channel_final]=optimum_model_function_2(a_store(secs(I_dr),section),b_store(secs(I_dr),section),alpha_store(secs(I_dr),section),Fs,35,channel_fil,0);
%               channel_final_35(section,:)=channel_final;
%             %  [channel_final_opt,aop,bop,alphaop]=optimum_model_function_3(channel_final,Fs);
%               if (length(channel_final)>1)  [T35_est_opt(section) EDT_35_est_opt(section) C80_35_est_opt(section) C50_35_est_opt(section) centre_35_est_opt(section) D_35_est_opt(section) ]=Room_acoustic_params_centre_ldr(channel_final,Fs,35);
%             end
%             else T35_est(section)=NaN;start_end_35(1:2)=NaN;
%               C80_35_est(section)=NaN;
%               C80_35_estm(section)=NaN;
%               C80_35_estmax(section)=NaN;
%               T35_est_opt(section)=nan;
%               EDT_35_est_opt(section)=nan;
%             end
%
%              I_dr=find(T25_est_sec>0);if size(T25_est_sec(I_dr)>0) [T25_est(section) I1]=min(T25_est_sec(I_dr));%start_end_25(1:2)=start_end(:,secs(I_dr(I1)));%%T25(section)=T25rt(I1);
%               C80_25_est(section)=min(C80_25(I_dr));C80_25_estm(section)=median(C80_25(I_dr));C80_25_estmax(section)=max(C80_25(I_dr));
%             centre_25_est(section)=min(centre_25(I_dr));centre_25_estm(section)=median(centre_25(I_dr));
%              % [channel_final]=optimum_model_function_orig(a_store(secs(I_dr)),b_store(secs(I_dr)),alpha_store(secs(I_dr)),Fs,25,channel_fil);
%             [channel_final]=optimum_model_function_2(a_store(secs(I_dr),section),b_store(secs(I_dr),section),alpha_store(secs(I_dr),section),Fs,25,channel_fil,0);
%             %lms_est=w_store(:)/max(abs(w_store(:)));
%             channel_fil=channel_fil/max(abs(channel_fil));
%             %Room_acoustic_params_plot(channel_final,Fs);
%             channel_final_25(section,:)=channel_final;
%             channel_final=channel_final;%(1:length(channel_fil));
%               [T25_est_opt(section) EDT_25_est_opt(section) C80_25_est_opt(section) C50_25_est_opt(section) centre_25_est_opt(section) D_25_est_opt(section) ]=Room_acoustic_params_centre_ldr(channel_final,Fs,25);
%              else T25_est(section)=NaN;start_end_25(1:2)=NaN;
%                  C80_25_est(section)=NaN;
%               C80_25_estm(section)=NaN;
%               C80_25_estmax(section)=NaN;
%              end
%
%              I_dr=find(T20_est_sec>0);if size(T20_est_sec(I_dr)>0) [T20_est(section) I1]=min(T20_est_sec(I_dr));%start_end_20(1:2)=start_end(:,secs(I_dr(I1)));%%T20(section)=T20rt(I1);
%                C80_20_est(section)=min(C80_20(I_dr));C80_20_estm(section)=median(C80_20(I_dr));C80_20_estmax(section)=max(C80_20(I_dr));
%                    centre_20_est(section)=min(centre_20(I_dr));centre_20_estm(section)=median(centre_20(I_dr));
%           [channel_final]=optimum_model_function_2(a_store(secs(I_dr),section),b_store(secs(I_dr),section),alpha_store(secs(I_dr),section),Fs,20,channel_fil,0);
%            %Room_acoustic_params_plot(channel_final,Fs);
%             channel_final_20(section,:)=channel_final;
%               [T20_est_opt(section) EDT_20_est_opt(section) C80_20_est_opt(section) C50_20_est_opt(section) centre_20_est_opt(section) D_20_est_opt(section) ]=Room_acoustic_params_centre_ldr(channel_final,Fs,20);
%
%
%              else T20_est(section)=NaN;start_end_20(1:2)=NaN;
%                  C80_20_est(section)=NaN;
%               C80_20_estm(section)=NaN;
%               C80_20_estmax(section)=NaN;
%              end
%
%             I_dr=find(T25_est_sec>0);if size(T25_est_sec(I_dr)>0)
%              T25_est_opt
%             channel_opt_store(section,1:length(channel_final_25))=channel_final_25(section,:);
%             else
%             channel_opt_store(section,1)=0;
%             end
% %filename=sprintf('%sMLE_3_track_%g_octave%g_all_picks_wind_%g_imp_.mat',name,section,Fc,imp)
% %filename=sprintf('MORE_PICKS/MLE_3_REAL_%s_track_%g_octave%g_all_picks_wind%g_maxRT_%g_imp_%g.mat',tag,track,Fc,wind,maxRT,imp)
% %save(filename,'channel_opt_store','a_store','b_store','alpha_store','len_store','DR','start_end','start_end2')
%
%
% end%%%End of section

% if Sections>1
%     II=find(channel_opt_store(:,1)>0)
%     if length(II)>1
%     channel_est=median(channel_opt_store(II,:));
%     else
%         channel_est=squeeze(channel_opt_store(II,:));
%     end
% else
%     channel_est=channel_opt_store;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This function is used to extracted the exponentional 
%      damping continuous segments of the speech signal  
%      Input:
%      x_E         ------   The envelope of the speech signal
%      Fs          ------   The sampling frequency
%      Best_s_num  ------   The numbers of the continuous segment which we want to extract
%
%      Output:
%      start_end   ------   The matrix which contains the start and end
%                           time index of the segments
%      An example of start_end: If we assume that Best_s_num = 3, which
%      means we want to extracted 3 continuous damping segments of the speech signal, then
%      start_end is a 2*3 matrix, where start_end(1,i) is the start time index of
%      the ith segments and start_end(2,i) is the end time index of the ith
%      segments
%      For example, start_end(1,1) = 2000, start_end(2,1) = 5000. Then the
%      first continuous damping segment is start from 2000 to 5000 (samples).
%                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function start_end = Choose_signal(x_E,Fs,Best_s_num,frame_t,to,step_size_s)
n = length(x_E);
%frame_t = 0.05;          % The length of the segment (seconds) which we want 
                         % to test weather it is a damping segments
                                
frame = floor(frame_t*Fs);  % The length of the segment (samples) which we want 
                            % to test weather it is a damping segments
                            
step_size = floor(step_size_s*Fs); % The step size when we choosing segments from the speech signal
                            % For example, if we choose one segment
                            % 0s-0.05s, then the next segment we want to
                            % test is then 0.01s-0.06s (step_size=0.01s)
% Hi Paul, you can modify frame_t and step_size according to your
% requirement
frame_num = floor((n-frame)/step_size+1);  % The number of segments which we want to test
kk = 1:frame;                             % The parameter which used in polyfit function

Var = zeros(frame_num,1);                 % This parameter is used to record the ratio of the polyfit line
                                          % If it is less than zero, then the segment envelope is a damping segment. 
y_tem = zeros(frame_num,2);               % This is used to record the result of polyfit
x_frame = zeros(frame,1);                 % This is used to record one segment
max_ytem = -log10(exp(6.91/to/Fs));        % This is the maximum value of the ratio of the polyfit line, with the RT = 3s
min_ytem = -log10(exp(6.91/0.005/Fs));    % This is the minimum value of the ratio of the polyfit line, with the RT = 0.005s
                                          % I haven't used them. 
min_ytem = -log10(exp(6.91/0.005/Fs));    
%6.91/(Fs*log(10^(-min_ytem)))=rt; 
max_count   = frame_num;
% h           = progressbar( [],0,'Polyfit progress' );
%  h = waitbar(0,'Please wait...performing Envelope segmentation');
count=0;
parfor i=1:frame_num 
%     count=count+1;
%     count/frame_num
%      (i/frame_num)
%     h = waitbar( h,1/frame_num );
%i/frame_num
    x_frame = x_E(((i-1)*step_size+1):((i-1)*step_size+frame) ); % Obtian one speech segment
     p= polyfit(kk',x_frame,1);                         % Using Polyfit function. you can check the help of polyfit
 y_tem(i,:)=p(:);
    store(i)=6.91/(Fs*log(10^(-p(1))));
    
    
    if((p<min_ytem) | (p>max_ytem))
  store(i)=0;
%    if(y_tem(i,1)>0)           % The ratio of the polyfit line. If it is larger than zero, then Var(i) = inf, and we will discard this segment later
        Var(i) = inf;
%         continue;
  end
    
end
% close(h)
%figure
%plot((0:(length(store)-1))/Fs,store)

[p_index] = Choose_seg(Var, step_size);      % Extracting damping segments according to Var, which records the ratios of the polyfit line
[pp,qq] = size(p_index);

out_num = min(pp,Best_s_num);
%out = zeros(length(x_E),out_num);
start_end = zeros(2,out_num);
for(i=1:out_num)
 
    [Best_continue,best_start_index] = max(p_index(:,2));  % Choosing the longest continuous damping segment
    best_start = p_index(best_start_index,1);
    p_index(best_start_index,2) = 0;                       % Remove the longest continuous damping segment, since it has been chosen
    start_end(1,i)= (best_start-1)*step_size+1;
    start_end(2,i)= (best_start-1)*step_size+frame+Best_continue*step_size;
end
%progressbar( h,-1 );





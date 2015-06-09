function [p_index] = Choose_p(p,step_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   p            records the ratio of the polyfit line of all segments
%
%   p_index      records the segment NO of the extracted continuous damping segments
%                For example if we have 50 segments, the NO is 1 to 50

%                p_index(i,1) denotes the NO of the first damping segments
%                             of ith continuous damping segments
%                p_index(i,2) denotes how many damping segments following 


%   An example of the result of p_index can be as follows:
%                p_index(1,1) = 2, p_index(1,2) = 6
%                p_index(2,1) = 10, p_index(2,2) = 16
%                p_index(3,1) = 19, p_index(3,2) = 20  ......
%
%                means the first continuous damping segment contains semgnt 2 to 6.
%                the second continuous damping segment contains semgnt 10 to
%                16.
%                the third continuous damping segment contains semgnt 19 to
%                20.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(p);
p_index=[];
k = 1;
pp = 1;
i=1;
while(i<n)
    k=1;
    if(p(i)<inf)  % p(i) records the ratio of the polyfit line for the ith segment
                  % if (p(i)<inf), then the ith segment is a damping segment
        p_index(pp,1) = i; % recording the number of the segments
        p_index(pp,2) = 0;
       
        while(((i+k)<=n) & (p(i+k)<inf) )
            p_index(pp,2) = p_index(pp,2)+1; % recording how many damping segments following
            k=k+1;
        end
        pp=pp+1;
        k = k+step_size;
      end
    i = i+k;
end
        
        



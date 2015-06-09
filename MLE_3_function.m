function [a,b,alpha]=MLE_3_function(x,Fs)
%%  Maximum likelihood fit of model to data

%define realistic poarameter search space using RTs
fromrt=(-6.91./log(0.99))/3000;
tort=(-6.91./log(0.99999))/3000;
from=exp(-6.91/(Fs*fromrt));
to=exp(-6.91/(Fs*tort));
steps1=5; %number of steps, more = more accurate
step1=(to-from)/(steps1-1);
expval=from:step1:to; %coarse search grid

%set up decay phase
y_re(1:length(x))=x./max(abs(x));
y_win=y_re;

% coarse search of parameter space
for bs=1:steps1
    b=expval(bs);
    ta=1;
    bb(bs)=b;
    a=from;
for as=1:steps1
    a=expval(as);
    aa(as)=a;
    alpha=0.5;
    
    %optimse wrt alpha parameter (i.e. convex weight of exponentials)
    [t_out1 fval]=fmincon(@(alpha) myfun_alpha(alpha,y_win,a,b),alpha,[],[],[],[],[0],[1],[],optimset('Algorithm' , 'active-set','LargeScale','off','Display','off','MaxFunEvals',300));
     alpha=t_out1;
    %sigma=sqrt(-sum(  (-1./(alpha*a.^i + (1-alpha)*b.^i).^2)   .*  (y_win.^2))/N);
    %L = - sum(log(alpha*a.^i + (1-alpha)*b.^i)) - sum( (   (alpha*a.^i + (1-alpha)*b.^i).^-2.*(y_win.^2) ) ./ (2*sigma^2) ) - N*log(2*pi*sigma^2)/2 ;%function
    test(as,bs)=-fval;
    test_alpha(as,bs)=alpha;

end

end
[vals,pos]=gmax(test);
a=aa(pos(1));
b=bb(pos(2));
alpha=test_alpha(pos(1),pos(2));

% fine search of parameter space, using a starting point from coarse search
x0=[a, b, alpha];
[x0 fval] = fmincon(@(x0) myfun_3(x0,y_win),x0,[],[],[],[],[from from 0],[to to 1],[],optimset('Algorithm' , 'active-set','LargeScale','off','MaxFunEvals',300,'Display','off'));
a=x0(1);
b=x0(2);
alpha=x0(3);
    
% figure
% semilogy(abs(y_win))
% hold on
% plot(alpha*a.^(1:length(y_win))+(1-alpha)*b.^(1:length(y_win)),'red')
% ylim([0.0001, 1])
function F = myfun_3(input,y_win)
% compute the log likelihood given the data in y_win
% this is for optimising wrt all parameters

a=input(1);
b=input(2);
alpha=input(3);

N=length(y_win);
i=((1:N)-1);
env=alpha*a.^i + (1-alpha)*b.^i;
sigma=sqrt(-sum(  (-1./(env).^2)   .*  (y_win.^2))/N);
L   = - sum(log(env)) - sum( (   (env).^-2.*(y_win.^2) ) ./ (2*sigma^2) ) - N*log(2*pi*sigma^2)/2 ;%function
F=-L;




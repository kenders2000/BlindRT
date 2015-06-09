function F = myfun(input,y_win,a,b)
% compute the log likelihood given the data in y_win
% this is for optimising wrt just alpha

alpha=input(1);
N=length(y_win);
i=((1:N)-1);
sigma=sqrt(-sum(  (-1./(alpha*a.^i + (1-alpha)*b.^i).^2)   .*  (y_win.^2))/N);
L       = - sum(log(alpha*a.^i + (1-alpha)*b.^i)) - sum( (   (alpha*a.^i + (1-alpha)*b.^i).^-2.*(y_win.^2) ) ./ (2*sigma^2) ) - N*log(2*pi*sigma^2)/2 ;%function
F=-L;







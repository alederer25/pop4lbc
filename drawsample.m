function [x,xi,xr]=drawsample(rsafe,Xss,gp,xref,tmax,bs)
% draw a new random scenario
% In:
%   rsafe      handle  performance specification
%   Xss        D  x 2  boundary of the considered state space
%   gp         object  GP class object
%   xref       handle  reference trajectory
%   tmax       1  x 1  final time for reference trajectory
%   bs         1  x 1  number of sampled scenarios
% Out:
%   x          D  x bs sampled states
%   xi         Dn x bs sampled uncertainties
%   xr         D  x bs sampled reference states
% Last edited: Armin Lederer, 04/2020

t=rand(1,bs)*tmax;
xr=xref(t);
x=xr; %initialization
for i=1:bs
    while(rsafe(x(:,i)-xr(:,i)))
        x(:,i)=rand(size(Xss,1),1).*(Xss(:,2)-Xss(:,1))+Xss(:,1);
    end
    [mu,sig]=gp.predict(x(:,i)');
    xi(:,i)=sig.*randn(length(mu),1);
end



end
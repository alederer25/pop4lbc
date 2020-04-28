function [kopt,X,Xi,Xr,isactive]=sop(epsbar,beta,Vdot_nom,Vdot_xi,Xss,r,gp,xref,tmax,batchsize,k0,ac)
% scenario optimization of the controller parameters
% In:
%   epsbar     1  x 1  chance constraint probability
%   beta       1  x 1  confidence level
%   Vdot_nom   handle  nominal Lyapunov function derivative
%   Vdot_xi    handle  uncertain Lyapunov function derivative component
%   Xss        D  x 2  boundary of the considered state space
%   r          handle  performance specificiation
%   gp         handle  probabilistic model
%   xref       handle  reference trajectory
%   tmax       1  x 1  final time considered for scenario optimization
%   batchsize  1  x 1  number of scenarios added per iteration
%   k0         N  x 1  initial controller parameters
%   ac         handle  cost function for scenario optimization
% Out:
%   kopt       N  x 1  optimized controller parameters
%   X          D  x M  sampled states
%   Xi         Dn x M  sampled model uncertainties
%   Xr         D  x M  sampled reference trajectory points
%   isactive   M  x 1  binary vector indicating the active constraints
% Last edited: Armin Lederer, 04/2020

epsilon=@(M,m)min(1,1-nthroot(beta/(M*nchoosek(M,m)),M-m));

M=0;
m=0;
epsi=1;

while(epsi>epsbar)
    M=M+1;
    [X(:,(M-1)*batchsize+1:M*batchsize),Xi(:,(M-1)*batchsize+1:M*batchsize),Xr(:,(M-1)*batchsize+1:M*batchsize)]=drawsample(r,Xss,gp,xref,tmax,batchsize);
    
    [kopt,isactive]=greedyopt(Vdot_nom,Vdot_xi,X,Xi,Xr,k0,ac);
    m=sum(isactive);
    epsi=epsilon(M*batchsize,m);
end

end
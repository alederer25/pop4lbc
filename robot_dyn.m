function y = robot_dyn(t,x,xp,n,tau,H,C,G)
% evaluate robot dynamics
% In:
%   t          1  x 1  time
%   x          D  x 1  state
%   xp         D  x 1  consistent initial condition
%   n          1  x 1  number of degrees of freedom  
%   tau        handle  joint torques
%   H          handle  mass matrix
%   C          handle  coriolis matrix
%   G          handle  gravity vector
% Out:
%   y          D  x 1  left side of differential equation
% Last edited: Armin Lederer, 04/2020
    y=zeros(2*n,1);
    y(1:n)=x(n+1:2*n)-xp(1:n);
    y(n+1:2*n)=H(x(1:n))*xp(n+1:end)+C(x(n+1:2*n),x(1:n))*x(n+1:2*n)+G(x(1:n))-tau(t,xp(n+1:end),x(n+1:2*n),x(1:n));

end
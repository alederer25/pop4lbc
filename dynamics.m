function y=dynamics(t,x,p_ff,p_fb,dyn,x_ref,xd_ref,xdd_ref,k)
% evaluate differential equation describing the dynamics of the controlled
% system
% In:
%   t          1  x 1  time
%   x          D  x 1  state
%   p_ff       handle  feedforward control component
%   p_fb       handle  feedback control component
%   dyn        handle  dynamics of the uncontrolled system
%   x_ref      handle  reference trajectory
%   xd_ref     handle  first derivative of the reference trajectory
%   xdd_ref    handle  second derivative of the reference trajectory
%   k          N  x 1  control gains
% Out:
%   y          N  x 1  left side of differential equation
% Last edited: Armin Lederer, 04/2020

y=zeros(2,1);
y(1)=x(2);
y(2)=dyn(p_ff(x(1),x(2),xdd_ref(t),xd_ref(t),x_ref(t))+p_fb(x-[x_ref(t);xd_ref(t)],k),x(1),x(2));

end
function [c,ceq]=constraint(k,Vnom,Vxi,X,Xi,Xr)
% computes the inequality constraints for the scenario optimization
% In:
%   k          1  x N  control parameters
%   Vnom       handle  nominal Lyapunov function derivative
%   Vxi        handle  uncertain component of the Lyapunov function derivative
%   X          D  x M  state samples
%   Xi         Dn x M  uncertainty samples
%   Xr         D  x M  reference trajectory samples
% Out:
%   c          1  x M  scenario constraints
%   ceq        0  x 0  equality constraints
% Last edited: Armin Lederer, 04/2020

ceq=[];

for i=1:size(X,2)
    c(i)=Vnom(X(:,i)-Xr(:,i),X(:,i),k)+Vxi(X(:,i)-Xr(:,i),Xi(:,i));
end

end
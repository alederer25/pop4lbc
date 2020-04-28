function [kopt,isactive]=greedyopt(Vdot_nom,Vdot_xi,X,Xi,Xr,k0,ac)
% optimization with greedy addition of scenario constraints
% In:
%   Vdot_nom   handle  nominal Lyapunov function derivative
%   Vdot_xi    handle  uncertain Lyapunov function derivative component
%   X          D  x M  state samples
%   Xi         Dn x M  uncertainty samples
%   Xr         D  x M  reference trajectory samples
%   k0         N  x 1  initial value of controller parameters
%   ac         handle  cost for scenario optimization
% Out:
%   kopt       N  x 1  optimized controller parameters
%   isactive   M  x 1  binary vector indicating the active constraints
% Last edited: Armin Lederer, 04/2020

vio=constraint(k0,Vdot_nom,Vdot_xi,X,Xi,Xr);
[~,idx]=max(vio);
isactive=false(size(X,2),1);
% ac=@(k)k'*k;
opt=optimset('Display','off','MaxFunEvals',5e4,'MaxIt',1e4,'Algorithm','sqp');
% opt=optimset('MaxFunEvals',1e5,'MaxIt',1e4,'Algorithm','sqp');
kopt=k0;

vio(idx)=1;%ensure at least one optimization

while(vio(idx)>1e-5)
    isactive(idx)=true;
    kinit=fmincon(@(k)1,k0,[],[],[],[],[],[],@(k)constraint(k,Vdot_nom,Vdot_xi,X(:,isactive),Xi(:,isactive),Xr(:,isactive)),opt);
    [kopt,~,flag]=fmincon(ac,kinit,[],[],[],[],[],[],@(k)constraint(k,Vdot_nom,Vdot_xi,X(:,isactive),Xi(:,isactive),Xr(:,isactive)),opt);
    if(flag<=0)
        [kopt,~,flag]=fmincon(@(k)ac(k)+sum(max(0,10000*constraint(k,Vdot_nom,Vdot_xi,X(:,isactive),Xi(:,isactive),Xr(:,isactive)))),kinit,[],[],[],[],1e-10*ones(size(k0)),[],@(k)constraint(k,Vdot_nom,Vdot_xi,X(:,isactive),Xi(:,isactive),Xr(:,isactive)),opt);
        if(flag<=0)
            flag
            error('scenario optimization might not converge due to infeasible constraints');
        end
    end
    vio=constraint(kopt,Vdot_nom,Vdot_xi,X,Xi,Xr);
    [~,idx]=max(vio);
end

end
function K_RBF = RBFs(k,q,qd,n_RBFs,n_dof)
% evaluate differential equation describing the dynamics of the controlled
% system
% In:
%   k          M  x 1  RBF parameters
%   q          D  x 1  joint angles
%   qd         D  x 1  joint velocities
%   n_RBFs     1  x 1  number of RBFs
%   n_dof      1  x 1  number of degrees of freedom
% Out:
%   K_RBF      N  x 1  control gain matrix
% Last edited: Armin Lederer, 04/2020
K_RBF=zeros(n_dof,n_dof);
for i=1:n_RBFs
    A = k(n_dof*(2+2*n_dof)*(i-1)+1:n_dof*(2+2*n_dof)*(i-1)+n_dof);
    L = k(n_dof*(2+2*n_dof)*(i-1)+n_dof+1:n_dof*(2+2*n_dof)*(i-1)+2*n_dof);
    C = k(n_dof*(2+2*n_dof)*(i-1)+2*n_dof+1:n_dof*(2+2*n_dof)*i);
    for j=1:n_dof
        K_RBF(j,j)=K_RBF(j,j)+A(j)*exp(-1./L(j)^2.*...
            sum((C(2*n_dof*(j-1)+1:2*n_dof*j)-[q;qd]).^2,1));
    end
end

K_RBF = K_RBF.^2+0.1*eye(n_dof);% add a small constant identity matrix to
%ensure positive definiteness of K_RBF

end
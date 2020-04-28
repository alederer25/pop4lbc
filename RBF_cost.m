function y=RBF_cost(k,n_RBFs,n_dof)
% evaluate cost for radial basis function network tuning
% In:
%   k          M  x 1  RBF network parameters
%   n_RBFs     1  x 1  number of RBFs
%   n_dof      1  x 1  number of degrees of freedom
% Out:
%   y          1  x 1  cost for the RBF parameters
% Last edited: Armin Lederer, 04/2020

y=0;
for i=1:n_RBFs
    A=k(n_dof*(2+2*n_dof)*(i-1)+1:n_dof*(2+2*n_dof)*(i-1)+n_dof,1); %amplitudes
    L=k(n_dof*(2+2*n_dof)*(i-1)+n_dof+1:n_dof*(2+2*n_dof)*(i-1)+2*n_dof,1);%length scales
    C=k(n_dof*(2+2*n_dof)*(i-1)+2*n_dof+1:n_dof*(2+2*n_dof)*i,1);%centers
    y=y+norm(L)^2+norm(A)^2;
end

end
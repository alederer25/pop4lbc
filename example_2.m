%% parameter definitions
clear all; 
clc;
warning('off'); %disable warnings arising from inaccuracy of the binomial coefficient computation 

Xss=[-1.2,1.2;-1.2,1.2;-1.2,1.2;-1.2,1.2]; %state space limitations
rsafe=@(e)(norm(e)<=0.1); %performance specification
tmax=20;%simulation time
epsbar=0.1; %violation probability
conf=1e-9; %confidence
batchsize=500; %number of new samples added in each iteration
Ntr=4; %number of training samples per dimension
sn = 0.1; %observation noise standard deviation
n_dof = 2; % Number of links
n_RBFs = 10; % number of radial basis functions per nonzero entry

% reference trajectory
q_ref=@(t)[sin(t);cos(t)];
qd_ref=@(t)[cos(t);-sin(t)];
qdd_ref=@(t)[-sin(t);-cos(t)];
xref=@(t)[q_ref(t);qd_ref(t)];

%initial state for simulation
q0=[0;1];
dq0=[1;0];
ddq0=[0;0];
xp0=[dq0;ddq0];
x0=[q0;dq0];

%initial controller parameters
for i=1:n_RBFs
    k0(n_dof*(2+2*n_dof)*(i-1)+1:n_dof*(2+2*n_dof)*(i-1)+n_dof,1)=randn(n_dof,1).^2; %amplitudes
    k0(n_dof*(2+2*n_dof)*(i-1)+n_dof+1:n_dof*(2+2*n_dof)*(i-1)+2*n_dof,1)=randn(n_dof,1).^2;%length scales
    k0(n_dof*(2+2*n_dof)*(i-1)+2*n_dof+1:n_dof*(2+2*n_dof)*i,1)=randn(2*n_dof^2,1);%centers
end

k0 = [k0;k0];


%% generate training data

%Robot model
m1=1;
m2=1;
l1=1;
l2=1;
g=10;

alpha=m1*(l1/2)^2+m2*(l1^2+(l2/2)^2);
beta=m2*l1*l2/2;
delta=m2*(l2/2)^2;
H=@(q) [alpha+2*beta*cos(q(2)), delta+beta*cos(q(2));delta+beta*cos(q(2)), delta];
C=@(q,qd) [-beta*sin(q(2))*qd(2), -beta*sin(q(2))*(qd(1)+qd(2)); beta*sin(q(2))*qd(1), 0];
G=@(q) [(m1+m2)*g*l1/2*sin(q(1))+m2*g*l2/2*sin(q(1)+q(2)); m2*g*l2/2*sin(q(1)+q(2))];

%prior model
rnd = @() 1+0.1*2*(0.5-rand);

m1h=rnd()*m1;
m2h=rnd()*m2;
l1h=rnd()*l1;
l2h=rnd()*l2;
gh=10;

alpha=m1h*(l1h/2)^2+m2h*(l1h^2+(l2h/2)^2);
beta=m2h*l1h*l2h/2;
delta=m2h*(l2h/2)^2;

hatH= H;
hatC=@(q,qd) [-beta*sin(q(2))*qd(2), -beta*sin(q(2))*(qd(1)+qd(2)); beta*sin(q(2))*qd(1), 0];
hatG=@(q) [(m1+m2)*g*l1/2*sin(q(1))+m2*g*l2/2*sin(q(1)+q(2)); m2*g*l2/2*sin(q(1)+q(2))];


X=ndgridj(-1*ones(2*n_dof,1),ones(2*n_dof,1),Ntr*ones(2*n_dof,1)); %training inputs

%training data generation
y=zeros(n_dof,Ntr^(2*n_dof));
ynom=zeros(n_dof,Ntr^(2*n_dof));
for i=1:Ntr^2
    y(:,i)=C(X(1:2,i),X(3:4,i))*X(3:4,i)+G(X(1:2,i));
    ynom(:,i)=hatC(X(1:2,i),X(3:4,i))*X(3:4,i)+hatG(X(1:2,i))+sn*randn(n_dof,1);
end

%% GP training and controller definition
gp1=fitrgp(X',(y(1,:)-ynom(1,:))','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact','Standardize',true,...
    'ConstantSigma',true,'Sigma',sn);
gp2=fitrgp(X',(y(2,:)-ynom(2,:))','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact','Standardize',true,...
    'ConstantSigma',true,'Sigma',sn);
gpModel=GPSSM(gp1,gp2);


p_ff=@(q,qd,qdd,qdd_ref,qd_ref,q_ref)hatH(q)*qdd_ref+hatC(qd,q)*qd_ref+hatG(q)+gpModel.predict([q;qd]'); %feedforward controller
p_fb=@(x,e,k)-RBFs(k(1:n_dof*(2+2*n_dof)*n_RBFs),x(1:n_dof),x(n_dof+1:2*n_dof),n_RBFs,n_dof)*e(1:n_dof)-...
    RBFs(k(n_dof*(2+2*n_dof)*n_RBFs+1:2*n_dof*(2+2*n_dof)*n_RBFs),x(1:n_dof),x(n_dof+1:2*n_dof),n_RBFs,n_dof)*e(n_dof+1:2*n_dof);%feedback controller

%% simulate system with high gain controller
Kp0= 1*eye(n_dof); %high gains
Kd0= 1*eye(n_dof);

u0=@(t,qdd,qd,q) p_ff(q,qd,qdd,qdd_ref(t),qd_ref(t),q_ref(t))-[Kp0, Kd0]*[q-q_ref(t); qd-qd_ref(t)];

opts = odeset('AbsTol',1e-3);
[t0,traj0] = ode15i(@(t,x,xp) robot_dyn(t,x,xp,n_dof,u0,H,C,G), [0,tmax],x0,xp0,opts);

%% define Lyapuonv function

epsilon=1e-1;

A = @(Kp,Kd,q,qd) [-Kd+epsilon*hatH(q), epsilon/2*(-Kd'+hatC(q,qd)); epsilon/2*(-Kd+hatC(q,qd)'), -epsilon*Kp];
Vdot_nom=@(e,x,k) [e(n_dof+1:2*n_dof); e(1:n_dof)]'*A(RBFs(k(1:n_dof*(2+2*n_dof)*n_RBFs),x(1:n_dof),x(n_dof+1:2*n_dof),n_RBFs,n_dof),...
    RBFs(k(n_dof*(2+2*n_dof)*n_RBFs+1:2*n_dof*(2+2*n_dof)*n_RBFs),x(1:n_dof),x(n_dof+1:2*n_dof),n_RBFs,n_dof),x(1:n_dof),x(n_dof+1:2*n_dof))*[e(n_dof+1:2*n_dof); e(1:n_dof)];
Vdot_xi=@(e,xi) [e(n_dof+1:2*n_dof); epsilon*e(1:n_dof)]'*[xi; xi];

%% run scenario optimization

[kopt,X,Xi,Xr,isactive]=sop(epsbar,conf,Vdot_nom,Vdot_xi,Xss,rsafe,gpModel,xref,tmax,batchsize,k0,@(k)RBF_cost(k,n_RBFs,n_dof));

%% simulate system with scenario optimized controller parameters

u=@(t,qdd,qd,q) p_ff(q,qd,qdd,qdd_ref(t),qd_ref(t),q_ref(t))+p_fb([q;qd],[q;qd]-[q_ref(t);qd_ref(t)],kopt);
[t,traj] = ode15i(@(t,x,xp) robot_dyn(t,x,xp,n_dof,u,H,C,G), [0,tmax],x0,xp0,opts);


%% plot results
figure(); 
plot(t,sqrt(sum((traj'-xref(t')).^2)));
hold on;
plot(t0,sqrt(sum((traj0'-xref(t0')).^2)))
title('control errors');


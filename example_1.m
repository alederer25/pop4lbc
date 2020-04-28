%% parameter definitions
clear all; 
clc;
warning('off'); %disable warnings arising from inaccuracy of the binomial coefficient computation 

Xss=[-1.2,1.2;-1.2,1.2]; %state space limitations
rsafe=@(e) (sum(e)<=0.1*sqrt(2)); %performance specification
tmax=2*pi;%simulation time
epsbar=0.01; %violation probability
beta=1e-9; %confidence
batchsize=200; %number of new samples added in each iteration
Ntr=21; %number of training samples per dimension
sn = 0.1; %observation noise standard deviation

% reference trajectory
q_ref=@(t)sin(t);
qd_ref=@(t)cos(t);
qdd_ref=@(t)-sin(t);
x_ref=@(t)[q_ref(t);qd_ref(t)];

%initial state for simulation
q0=0;
qd0=1;



%% generate training data

c=rand(1,1)*2*pi; %random parameter of dynamics
dyn=@(u,q,qd)u-qd-1-(qd^2*sin(q-c)-sin(c))/(cos(q-c)-1.1/cos(q-c)); %unknown dynamics

X=ndgridj([-1;-1],[1;1],Ntr*ones(2,1)); %training inputs

model=@(q,qd,tau)tau-1-qd; %prior model

%training data generation
y=zeros(1,Ntr^2);
ynom=zeros(1,Ntr^2);
for i=1:Ntr^2
    y(i)=dyn(0,X(1,i),X(2,i))+randn(1)*sn;
    ynom(i)=model(X(1,i),X(2,i),0);
end

%% GP training and controller definition
gp=fitrgp(X',(y-ynom)','KernelFunction','ardsquaredexponential','FitMethod','exact','PredictMethod','exact','Standardize',true,...
    'ConstantSigma',true,'Sigma',sn);


p_ff=@(q,qd,qdd_ref,qd_ref,q_ref)1*qdd_ref-model(q,qd_ref,0)-1*gp.predict([q;qd]'); %feedforward controller
p_fb=@(e,k)-k*sum(e);%feedback controller

%% simulate system with high gain controller
k=10; %high gain controller gains

[t0,traj0]=ode45(@(t,x)dynamics(t,x,p_ff,p_fb,dyn,q_ref,qd_ref,qdd_ref,k),[0,2*pi],[q0;qd0]); %simulate high gain trajectory


%get control inputs along simulated trajectory
u0=zeros(1,size(t0,1));
for i=1:length(t0)
    u0(i)=p_ff(traj0(i,1),traj0(i,2),qdd_ref(t0(i)),qd_ref(t0(i)),q_ref(t0(i)))+p_fb(traj0(i,:)'-x_ref(t0(i)),k);
end

%% define Lyapuonv function

Vdot_nom=@(e,x,k)-k*sum(e)^2; %nominal Lyapunov derivative component
Vdot_xi=@(e,xi)sum(e)*xi; %uncertain Lyapunov derivative component


%% run scenario optimization
k0=1;
ac=@(k)k^2;
[kopt,X,Xi,Xr,isactive]=sop(epsbar,beta,Vdot_nom,Vdot_xi,Xss,rsafe,gp,x_ref,tmax,batchsize,k0,ac);

%% simulate system with scenario optimized controller gains
[t,traj]=ode45(@(t,x)dynamics(t,x,p_ff,p_fb,dyn,q_ref,qd_ref,qdd_ref,kopt),[0,tmax],[q0;qd0]);

u=zeros(1,size(t,1));
for i=1:length(t)
    u(i)=p_ff(traj(i,1),traj(i,2),qdd_ref(t(i)),qd_ref(t(i)),q_ref(t(i)))+p_fb(traj(i,:)'-x_ref(t(i)),k);
end


%% plot results
figure(); 
plot(t,sqrt(sum((traj'-x_ref(t')).^2)));
hold on;
plot(t0,sqrt(sum((traj0'-x_ref(t0')).^2)))
title('control errors');

figure(); 
plot(t,u);
hold on;
plot(t0,u0);
title('control inputs');

disp(['scenario optimized control gain: ', num2str(kopt(1))]);
disp(['fixed high gain: ', num2str(k(1))]);



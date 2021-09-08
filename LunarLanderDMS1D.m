clc;
clear

ge=9.8;
g=ge/6;
Isp=310;

L = @(x,u,t)(0);% lagrange performance index
M = @(x,T)(-x(3));% meyer Performance index
% state control contrastints
scon = @(x,u)[ u - 7500 ; 1500 - u; -x(1) ;-x(2)];
% terminal constraint
psi = @(x,T) [ x(1) ; x(2)];
% system ODE equation
f_ode = @(x,u,t)[x(2);
            -g+u/x(3);
            -u/(ge*Isp)];
x_0 = [ 200 ; -20 ; 1200];

tf = 30;
Nodes = 30;
m = 1; % number of control input

% g.t_guess=zeros(1,Nodes+1);
% g.x_guess=zeros(n,Nodes+1);
% g.u_guess=zeros(p,Nodes+1);
guess.t_guess=linspace(0,tf,Nodes+1);
guess.x_guess=[linspace(x_0(1),0,Nodes+1);
            linspace(x_0(2),0,Nodes+1);
            linspace(x_0(3),1000,Nodes+1)];
guess.u_guess=[guess.x_guess(3,:).*([diff(guess.x_guess(2,:)),0]+g)];

[X,U,t,J] = DMS(L,M,scon,psi,f_ode,x_0,m,tf,Nodes,guess); % or DMS
%[X,U,t,J] = DMS(L,M,h,r,f,x_0,m,T,N); % or DMS
X=X';t=t';U=U';
%%
figure(1)
subplot(3,1,1)
plot(t,X(:,1)),ylabel('height/m','Interpreter','latex');
subplot(3,1,2)
plot(t,X(:,2)),ylabel('velocity (m/s)','Interpreter','latex');
subplot(3,1,3)
plot(t,X(:,3)),ylabel('mass/kg','Interpreter','latex');
xlabel('time/s','Interpreter','latex');

figure(2)
plot(t,U),ylabel('thrust/N','Interpreter','latex');
xlabel('time/s','Interpreter','latex');

                                                                      clc;
clear

ge=9.8;
g=ge/6;
Isp=310;

L = @(x,u,t)(0);% lagrange performance index
M = @(x,T)(-x(3));% meyer Performance index
% state control contrastints
u_min = 1500;
u_max = 7500;
scon = @(x,u)[ u - u_max ; u_min - u ];
% terminal constraint
psi = @(x,T) [ x(1) ; x(2)];
% system ODE equation
f_ode = @(x,u,t)[x(2);
            -g+u/x(3);
            -u/(ge*Isp)];
x_0 = [ 200 ; -2 ; 1200];

tf = 30;
Nodes = 50;
m = 1; % number of control input

% g.t_guess=zeros(1,Nodes+1);
% g.x_guess=zeros(n,Nodes+1);
% g.u_guess=zeros(p,Nodes+1);
guess.t_guess=linspace(0,tf,Nodes+1);
vel_guess = linspace(x_0(2),-1,Nodes+1);
guess.x_guess=[linspace(x_0(1),0,Nodes+1);
            vel_guess;
            linspace(x_0(3),1150,Nodes+1)];
% guess.u_guess=[x_0(3).*([0,diff(guess.x_guess(2,:))]+g)];

guess.u_guess=[zeros(1,40),u_min*ones(1,Nodes-40+1)];

[X,U,t,J] = DSS(L,M,scon,psi,f_ode,x_0,m,tf,Nodes,guess); % or DMS
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
plot(t,U,'b',[0,tf]',[u_min,u_max;u_min,u_max],'k--'),
legend('u','u_{min}','u_{max}')
ylabel('thrust/N','Interpreter','latex');
xlabel('time/s','Interpreter','latex');

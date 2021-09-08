clc;
clear

L = @(x,u,t)(1);% lagrange performance index
M = @(x,T)(0);% meyer Performance index
% state control contrastints
scon = @(x,u)[  u - 1.5*pi ;
                0 - u];
            
% scon = @(x,u)[ -x(1);
%                 -x(2);
%                 -x(3);
%                 -x(4)-10;
%                 x(1)-30;
%                 x(2)-30;
%                 x(3)-20;
%                 x(4)-10;
%                 u - pi ;
%                 -pi - u];
% terminal constraint
psi = @(x,T) [ x(1)-5 ;
                x(2)-5;
                x(3);
                x(4)-5];
% system ODE equation
f_ode = @(x,u,t)[x(3);
                x(4);
            3*cos(u);
            3*sin(u)];
        
x_0 = [ 0;0;5;0];

tf = 1.6;
Nodes = 100;
m = 1; % number of control input

% g.t_guess=zeros(1,Nodes+1);
% g.x_guess=zeros(n,Nodes+1);
% g.u_guess=zeros(p,Nodes+1);
tg=linspace(0,tf,Nodes+1); guess.t_guess=tg;
radius=5;
omega=1;
guess.x_guess=[radius*sin(omega*tg);
            radius*(1-cos(omega*tg));
            radius*omega*cos(omega*tg);
            radius*omega*(sin(omega*tg))];
%         temp=[diff(guess.x_guess(2,:));diff(guess.x_guess(1,:))];
%         temp=[temp,temp(:,end)];
% guess.u_guess=atan2(temp(1,:),temp(2,:));
guess.u_guess=[ones(1,Nodes+1)*pi];

[X,alpha,t,J] = DMS(L,M,scon,psi,f_ode,x_0,m,[],Nodes,guess); % or DMS
%[X,U,t,J] = DMS(L,M,h,r,f,x_0,m,T,N); % or DMS
X=X';t=t';alpha=alpha';
%%
figure(1)
subplot(2,1,1),title('velocity');
plot(t,X(:,3)),ylabel('$v_x$ (m/s)','Interpreter','latex');
subplot(2,1,2)
plot(t,X(:,4)),ylabel('$v_y$ (m/s)','Interpreter','latex');

figure(2)
plot(X(:,1),X(:,2)),title('trajectory')
ylabel('y/m','Interpreter','latex');xlabel('x/m','Interpreter','latex');

figure(3)
plot(t,alpha),ylabel('$\theta$/rad','Interpreter','latex');
xlabel('time/s','Interpreter','latex');
title('turning angle \alpha');
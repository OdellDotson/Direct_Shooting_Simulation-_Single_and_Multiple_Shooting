clc;
clear
L = @(x,u,t)(u*u);
M = @(x,T)(0);
u_max = 10;
h = @(x,u)[ u - 10 ; -10 - u];
e = 0.0001;
r = @(x,T) [ x(1) ; x(2)];


m = 0.1;
l = 1;
g = 9.81;
b = 0.1;

f = @(x,u,t)[x(2);
    (m*g*l*sin(x(1))-b*x(2)+u)/(m*l^2)];
x_0 = [ pi/2 ; 0 ];
T = 3;
N = 30;

m = 1; % number of control input

[X,U,t,J] = DSS(L,M,h,r,f,x_0,m,T,N); % or DMS
% (L,M,scon,psi,f_ode,x_0,p,tf,Nodes,guess)
%[X,U,t,J] = DMS(L,M,h,r,f,x_0,m,T,N); % or DMS
%%

subplot(2,1,1)
plot(t,X),xlabel('time/s'),y;ane;('states');
subplot(2,1,2)
plot(t,U),xlabel('time/s'),y;ane;('control');

figure
hold on;
t_previous = -T/N;
for i = 1:size(X,2)
    
    p = plot([0 l*sin(X(1,i))],[0 l*cos(X(1,i))],'b','LineWidth',3);
    axis([-1 1 -1 1])
    axis equal
    
    pause(t(i)-t_previous);
    t_previous = t(i);
    if i ~= size(X,2)
        delete(p);
    end
end

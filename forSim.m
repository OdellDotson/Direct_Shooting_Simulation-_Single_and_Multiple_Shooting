function [X,ts] = forSim(f,x_0,U,T,N)
% u(x,t)

ts = linspace(0,T,N+1);
X = zeros(numel(x_0),N+1);
dt = 0.1*T/N;
X(:,1) = x_0;
for i = 1:N
    tnow=ts(i);
    Xtemp=X(:,i);
    while tnow < ts(i+1)
        Xtemp = rk4(f,Xtemp,U(:,i),tnow,dt);
        tnow=tnow+dt;
    end
    X(:,i+1) = Xtemp;
end

end


function x_next = rk4(f,x,u,t,dt)
k1 = f(x,u,t);
k2 = f(x+k1*dt/2,u,t+dt/2);
k3 = f(x+k2*dt/2,u,t+dt/2);
k4 = f(x+k3+dt/2,u,t+dt);

x_next = x + (k1+2*k2+2*k3+k4)*dt/6;


end

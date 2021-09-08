function [X,U,t,J] = DMS(L,M,scon,psi,f_ode,x_0,m,tf,Nodes,guess)

% initialize
n = numel(x_0);
if(nargin==9)
    w_0 = 0.01*ones((n+m)*(Nodes+1),1);
else
    t_guess=guess.t_guess;
    x_guess=guess.x_guess;% n*(nodes+1)
    u_guess=guess.u_guess;
    w_0=[reshape(x_guess,[1,n*(Nodes+1)]),...
         u_guess]';
end
% % w_0 = 2*rand(m*(N+1),1)-1;
if isempty(tf)
    tf0=1.6;
    w_0 = [w_0; tf0];
end

% defining objective function
fun = @(w)bolzaDMS(w,L,M,f_ode,x_0,m,tf,Nodes);

% defining constraints
nonlcon = @(w)constraintDMS(w,f_ode,scon,psi,x_0,m,tf,Nodes);

% solving transformed NLP
opts = optimoptions(@fmincon,'TolX',1e-3,'TolFun',1e-4,'GradObj','on',...
    'TolCon',1e-2,'GradConstr','off','MaxFunEvals',10000);
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon,opts);


% extract U and X
if isempty(tf)
    tf = w(end);
    w = w(1:end-1);
end
X = reshape(w(1:n*(Nodes+1)),n,Nodes+1);
U = reshape(w(n*(Nodes+1)+1:(n+m)*(Nodes+1)),m,Nodes+1);
U = [U(:,1:end-1) U(:,end-1)];
t = linspace(0,tf,Nodes+1);
end

function [J,pJpw] = bolzaDMS(w,L,M,f_ode,x_0,m,tf,Nodes)
    lw=length(w);
    % extract U
    if isempty(tf)
        lw=lw-1;
        tff = w(end);
        w = w(1:lw);
    end
    % extract U
    n = numel(x_0);
    X_guess = reshape(w(1:n*(Nodes+1)),n,Nodes+1);
    U_guess = reshape(w(n*(Nodes+1)+1:(n+m)*(Nodes+1)),m,Nodes+1);

    % numerically integrate the Lagrange term
    L_a = @(Lambda,v,t) augment(L,v,t,n);

    v = [X_guess;U_guess];

    dt = tff/Nodes;
    int_Lambda=[0,0];
    for i = 1:Nodes
        int_Lambda = forSim(L_a,int_Lambda(2),v(:,i),dt,1);
    end
    Lambda =  int_Lambda(end);
    % evaluate total Bolza objective cost
    J = Lambda + M(X_guess(:,end),tff);

    if nargout > 1 % gradient required
        pJpw = [zeros(lw,1);
                   1 ];
    end
    % augmented Lagrange term to work with forSim
        function L_a = augment(L,v,t,n)
            x = v(1:n);
            u = v(n+1:end);
            L_a = L(x,u,t);
        end

end

function [inCon,eqCon] = constraintDMS(w,f_ode,scon,psi,x_0,m,tf_in,Nodes)
    % inCon < 0
    lw=length(w);
    % extract U
    if isempty(tf_in)
        lw=lw-1;
        tf_dec = w(end);
        w = w(1:lw);
        isTFree = 1;
    end
    n = numel(x_0);
    X_try = reshape(w(1:n*(Nodes+1)),n,Nodes+1);
    U_guess = reshape(w(n*(Nodes+1)+1:(n+m)*(Nodes+1)),m,Nodes+1);

    % allocate space for inCOn
    % total number of discretized path constraints 
    n_h = size(scon(X_try(:,1),U_guess(:,1)),1); 
    inCon = zeros(n_h*(Nodes+1),1);
    % total number of discretized path constraints 
    n_g = size(X_try(:,1),1); 
    eqCon = zeros(n_g*(Nodes+1),1);
    dt = tf_dec/Nodes;

    eqCon(1:n,1) = X_try(:,1) - x_0;
    for i = 1:Nodes
        % state and control path constraints
        inCon(n_h*(i-1)+1:n_h*i) = scon(X_try(:,end),U_guess(:,i));
        X_sim = forSim(f_ode,X_try(:,i),[U_guess(:,i) U_guess(:,i)],dt,1);

        % equality constraints across every phases
        eqCon(i*n+1:(i+1)*n,1)=X_try(:,i+1)-X_sim(:,2);
    end
    eqCon(end+1:end+size(psi(X_try(:,end),tf_dec),1)) = psi(X_try(:,end),tf_dec);

    % fill in this part
    %inCon(end+1:end+size(r(X_sim(:,end),T),1)) = ;

    if isTFree
        inCon(end+1) = -tf_dec + 0;% tf > 0
    end

end
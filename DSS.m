function [X,U,t,J] = DSS(L,M,scon,psi,f_ode,x_0,m,tf,Nodes,guess)
% f
% initialize
n = numel(x_0);
if(nargin==9)
    w_0 = ones((m)*(Nodes+1),1);
else
    t_guess=guess.t_guess;
    x_guess=guess.x_guess;% n*(nodes+1)
    u_guess=guess.u_guess;
    w_0=[u_guess]';
end
if isempty(tf)
    tf0=2;
    w_0 = [w_0; tf0];
end

% defining objective function
fun = @(w)bolzaDSS(w,L,M,f_ode,x_0,m,tf,Nodes);

% defining constraints
nonlcon = @(w)constraintDSS(w,f_ode,scon,psi,x_0,m,tf,Nodes);

% solving transformed NLP
opts = optimoptions(@fmincon,'Display','iter','MaxFunEvals',10000,...
    'TolX',1e-2,'TolFun',1e-4,'GradObj','off',...
    'TolCon',1e-3,'GradConstr','off');
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon,opts);

% extract U and X
if isempty(tf)
    tf = w(end);
    w = w(1:end-1);
end

U = reshape(w,m,numel(w)/m);
U = [U(:,1:end-1) U(:,end-1)];


% fill in this part
[X,t] = forSim(f_ode, x_0, U, tf, Nodes); %discretized state trajectory, and time
end

function [J,pJpw] = bolzaDSS(w,L,M,f_ode,x_0,m,tf,Nodes)
    lw=length(w);
    % extract U
    if isempty(tf)
        lw=lw-1;
        tff = w(end);
        w = w(1:lw);
    else
        tff=tf;
    end
    U_guess = reshape(w,m,lw/m);

    [X_sim,~] = forSim(f_ode, x_0,U_guess, tff, Nodes);

    % numerically integrate the Lagrange term
    n = numel(x_0);
    L_a = @(Lambda,v,t) augment(L,v,t,n); % augmented Lagrange term
    v = [X_sim;U_guess];

    Lambda = forSim(L_a,0,v,tff,Nodes);
    Lambda = Lambda(end); % integrate Lagrange term

    % evaluate total Bolza objective cost
    J = Lambda + M(X_sim(:,end),tff);

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

function [inCon,eqCon,gradc,gradceq] = constraintDSS(w,f_ode,scon,psi,x_0,m,tf_in,Nodes)
    % inCon < 0
    lw=length(w);
    % extract U
    if isempty(tf_in)
        lw=lw-1;
        tf_dec = w(end);
        w = w(1:lw);
        isTFree = 1;
    else
        tf_dec=tf_in;
        isTFree = 0;
    end
    U_guess = reshape(w,m,lw/m);

    % forward simulate to get X
    [X_sim,~] = forSim(f_ode, x_0,U_guess, tf_dec, Nodes);

    % allocate space for inCOn
    n_h = size(scon(X_sim(:,1),U_guess(:,1)),1); % total number of discretized path constraints
    inCon = zeros(n_h*(Nodes+1),1);

    % fill in the inequality constrants
    for i = 1:Nodes+1
        % fill in this part
        inCon(n_h*(i-1)+1:n_h*i) = scon(X_sim(:,end),U_guess(:,i));
    end

    %inCon(end+1:end+size(r(X_sim(:,end),T),1)) = r;

    if isTFree
        inCon(end+1) = -tf_dec + 0;
    end

    % terminal constraints
    eqCon = psi(X_sim(:,end),tf_dec);
    
    % gradient of equality or inequality constraints
    if nargout>2
        gradc
        gradceq
    end
end
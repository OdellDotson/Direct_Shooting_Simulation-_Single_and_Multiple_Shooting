function [X,U,t,J] = DSS(L,M,h,r,f,x_0,m,T,N)

% initialize
w_0 = 0.01*ones(m*(N+1),1);
if isempty(T)
    w_0 = [w_0; 10];
end

% defining objective function
fun = @(w)bolzaDSS(w,L,M,f,x_0,m,T,N);

% defining constraints
nonlcon = @(w)constraintDSS(w,f,h,r,x_0,m,T,N);

% solving transformed NLP
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon);

% extract U and X
if isempty(T)
    T = w(end);
    w = w(1:end-1);
end

U = reshape(w,m,numel(w)/m);
U = [U(:,1:end-1) U(:,end-1)];


% fill in this part
[X,t] = forSim(f, x_0,U, T, N); %discretized state trajectory, and time

end

function J = bolzaDSS(w,L,M,f,x_0,m,T,N)

% extract U
if isempty(T)
    T = w(end);
    w = w(1:end-1);
end
U_guess = reshape(w,m,numel(w)/m);

% fill in this part
[X_sim,~] = forSim(f, x_0,U_guess, T, N);

% numerically integrate the Lagrange term
n = numel(x_0);
L_a = @(Lambda,v,t) augment(L,v,t,n); % augmented Lagrange term
v = [X_sim;U_guess];

Lambda = forSim(L_a,0,v,T,N); 
Lambda = Lambda(end); % integrate Lagrange term

% evaluate total Bolza objective cost

% fill in this part
J = Lambda + M(X_sim(:,end),T);

% augmented Lagrange term to work with forSim
    function L_a = augment(L,v,t,n)
        x = v(1:n);
        u = v(n+1:end);
        L_a = L(x,u,t);
    end

end

function [inCon,eqCon] = constraintDSS(w,f,h,r,x_0,m,T,N)

% extract U
isTFree = 0;
if isempty(T)
    T = w(end);
    w = w(1:end-1);
    isTFree = 1;
end
U_guess = reshape(w,m,numel(w)/m);

% forward simulate to get X
% Fill in this part
[X_sim,~] = forSim(f, x_0,U_guess, T, N);

% allocate space for inCOn
n_h = size(h(X_sim(:,1),U_guess(:,1)),1); % total number of discretized path constraints 
inCon = zeros(n_h*(N+1),1);

% fill in the inequality constrants
for i = 1:N+1,
    % fill in this part
    inCon(n_h*(i-1)+1:n_h*i) = h(X_sim(:,end),U_guess(:,i));
end

% fill in this part
%inCon(end+1:end+size(r(X_sim(:,end),T),1)) = r;

if isTFree
    inCon(end+1) = -T;
end

% terminal constraints
% fill in this part
eqCon = r(X_sim(:,end),T);

end
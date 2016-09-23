function [X,U,t,J] = DMS(L,M,h,r,f,x_0,m,T,N)

% initialize
n = numel(x_0);
w_0 = 0.01*ones((n+m)*(N+1),1);
% w_0 = 2*rand(m*(N+1),1)-1;
if isempty(T)
    w_0 = [w_0; 10];
end

% defining objective function
fun = @(w)bolzaDMS(w,L,M,f,x_0,m,T,N);

% defining constraints
nonlcon = @(w)constraintDMS(w,f,h,r,x_0,m,T,N);

% solving transformed NLP
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon);

% extract U and X
if isempty(T)
    T = w(end);
    w = w(1:end-1);
end
X = reshape(w(1:n*(N+1)),n,N+1);
U = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);
U = [U(:,1:end-1) U(:,end-1)];
t = linspace(0,T,N+1);
end

function J = bolzaDMS(w,L,M,f,x_0,m,T,N)

% extract U
n = numel(x_0);
X_guess = reshape(w(1:n*(N+1)),n,N+1);
U_guess = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);

if isempty(T)
    T = w(end);
end
% numerically integrate the Lagrange term
L_a = @(Lambda,v,t) augment(L,v,t,n);

v = [X_guess;U_guess];

dt = T/N;
Lambda = 0;
for i = 1:N
    int_Lambda = forSim(L_a,Lambda,v(:,i),dt,1);
    Lambda = Lambda + int_Lambda(end);
end

% evaluate total Bolza objective cost
%Fill in this part
J = Lambda + M(X_guess(:,end),T);

% augmented Lagrange term to work with forSim
    function L_a = augment(L,v,t,n)
        x = v(1:n);
        u = v(n+1:end);
        L_a = L(x,u,t);
    end

end

function [inCon,eqCon] = constraintDMS(w,f,h,r,x_0,m,T,N)

% extract U
isTFree = 0;
if isempty(T)
    T = w(end);
    isTFree = 1;
end
n = numel(x_0);
X_guess = reshape(w(1:n*(N+1)),n,N+1);
U_guess = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);

% allocate space for inCOn
n_h = size(h(X_guess(:,1),U_guess(:,1)),1); % total number of discretized path constraints 
inCon = zeros(n_h*(N+1),1);

n_g = size(X_guess(:,1),1); % total number of discretized path constraints 
eqCon = zeros(n_g*(N+1),1);
dt = T/N;
% fill in the inequality constrants
for i = 1:N,
    % fill in this part
    inCon(n_h*(i-1)+1:n_h*i) = h(X_guess(:,end),U_guess(:,i));
    X_sim = forSim(f,X_guess(:,i),[U_guess(:,i) U_guess(:,i)],dt,1);
    % fill in this part
    eqCon(end+1:end+size(r(X_guess(:,end),T),1)) = r(X_guess(:,end),T);
    
    eqCon(end+1:end+n) = X_guess(:,1) - x_0;
end

% fill in this part
%inCon(end+1:end+size(r(X_sim(:,end),T),1)) = ;

if isTFree
    inCon(end+1) = -T;
end

end
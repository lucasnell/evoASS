function [V,N,U,P,WN,WP] = doManyTimeSteps(V,N,U,P,T,q,mutProb,SDmut)

% -----------------
% Adaptive dynamics version
% -----------------


etaN = 0.1;
etaP = 0.1;

r = .1;
a = 1;
f = 0.1;
cc = 1;
m = .01;
g = .01;
b = .05;% 0.05 or 0.5
minpopsize = 1e-4;

C = etaN*ones(q) + (1-etaN)*eye(q);
D = etaP*ones(q) + (1-etaP)*eye(q);

A = Afunc(V,C,a,f);% specific density dependence, increases with V
B = Bfunc(V,U,b); % attack rates, declines with V^2 and U^2. B(i,j) is for prey i, pred j
M = Mfunc(U,D,m,g); % predator mortality, decreases with U

for t=1:T
    % first one ecological time step:
	WN = r*(1 - A*sum(N) - B * P);
    WP = cc*B'*N - M;
    Nt = N.*exp(WN);
    Pt = P.*exp(WP);
    N = Nt;
    P = Pt;
    % extinctions:
    Nextinct = N < minpopsize;
    if any(Nextinct)
        V(Nextinct,:) = [];
        A(Nextinct) = [];
        B(Nextinct,:) = [];
        N(Nextinct) = [];
        WN(Nextinct) = [];
        if isempty(N)
            disp('N extinction!')
            break
        end
    end
    Pextinct = P < minpopsize;
    if any(Pextinct)
        U(Pextinct,:) = [];
        B(:,Pextinct) = [];
        M(Pextinct) = [];
        P(Pextinct) = [];
        WP(Pextinct) = [];

        if isempty(P)
            disp('P extinction!')
            break
        end
    end
    % prey mutations:
    Nmut = rand(size(N)) < mutProb*N;
    if any(Nmut)
        for i=find(Nmut)'
            V(end+1,:) = max(0, V(i,:) + SDmut*randn(1,q) );
            N(end+1,1) = 1.01*minpopsize;
            A =  Afunc(V,C,a,f); % specific density dependence, increases with V
            B = Bfunc(V,U,b); % attack rates, declines with V^2 and U^2. B(i,j) is for prey i, pred j
        end
    end
    % pred mutations:
    Pmut = rand(size(P)) < mutProb*P;
    if any(Pmut)
        for i=find(Pmut)'
            U(end+1,:) = max(0, U(i,:) + SDmut*randn(1,q) );
            P(end+1,1) = 1.01*minpopsize;
            B = Bfunc(V,U,b); % attack rates, declines with V^2 and U^2. B(i,j) is for prey i, pred j
            M = Mfunc(U,D,m,g); % predator mortality, decreases with U
        end
    end
end


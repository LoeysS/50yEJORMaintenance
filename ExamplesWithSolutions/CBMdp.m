function [ M , g , P0, P1] = CBMdp(pmf,Cu,Cp,epsilon)

L = length(pmf);
% compute transition probabilities
q = pmf;

P0 = zeros(L,L);
for i=1:L
    P0(i,i:end-1) = q(1:L-i);
    P0(i,end) = 1 - sum(P0(i,1:end-1));
end
P1 = zeros(L,L);
for i=1:L
    P1(i,:) = q(1:L);
    P1(i,end) = 1 - sum(q(1:L-1));
end
V1 = zeros(L,1); %states and stages to go
V = zeros(L,1);
V(L) = Cu;
term1 = 0; % place holders for minimization terms
term2 = 0;
%M = L*ones(1,N+1);
span = 2*epsilon;
while span >= epsilon
    % go through state space
    for x=1:L-1
        term1 = Cp + P1(x,:)*V;
        term2 = P0(x,:)*V;
        V1(x) = min(term1,term2);
    end
    V1(L) = Cu + P1(L,:)*V;
    span = max(V1-V) - min(V1-V);
    g = ( max(V1-V) + min(V1-V) ) / 2;
    V = V1;
end

% determine threshold
M=L-1;
for x=1:L
    if Cp + dot(P1(x,:),V) <= dot(P0(x,:),V)
        M = x-1;
        break
    end
end

end



function [CostRate, T] = cbmBayesPoisson(Cu,Cp,L,a,b)

%create three dimensional array with transition probabilities

% compute appropriate level to truncate state space in t direction
tmax = ceil(2*L/gaminv(0.001,a,1/b));   
xmax = L;
% create vectors of the possible r and p parameters of the negative
% binomial distribution that may occur as you go through the state space
rvec = (a+(0:xmax))';
pvec = ( (b + (0:tmax)) ./ (b + (0:tmax) + 1.0))';
% create matrices with relevant parameters to compute transition
% probabilities in one call (vectorize).
xarg = zeros([xmax+1,L+1,tmax+1]);
rarg = zeros([xmax+1,L+1,tmax+1]);
parg = zeros([xmax+1,L+1,tmax+1]);
for x=0:L
    for ip=1:length(pvec)
        for ir=1:length(rvec)
            xarg(x+1,ir,ip)=x;
            parg(x+1,ir,ip)=pvec(ip);
            rarg(x+1,ir,ip)=rvec(ir);
        end
    end
end
% compute all possible transition probabilities.
pmf = nbinpdf(xarg,rarg,parg);
% normalize transition probabilities after truncation
pmfnew = pmf(:,1,1);
pmfnew(end) = 1.0 - sum(pmfnew(1:end-1));

% allocate memory for value functions
Vn = zeros(xmax+1,tmax);
Vn1 = zeros(xmax+1,tmax);

span = 1.1;
T = zeros(tmax,1);

iter = 0;
% loop over all states
while span > 1e-6
    iter = iter + 1;
    for t=1:tmax
        T(t) = 0.0;        
        for x = 0:xmax

            if( x==xmax )
                Vn1(x+1,t) = Cu + dot(pmfnew,Vn(:,1));
            else
                %cost of preventive replacement
                replace = Cp + dot(pmfnew,Vn(:,1));
            end
            
            if ( x<xmax && t<tmax )
                %compute cost of no replacement
                notreplace = (pmf(1:(xmax-x),x+1,t+1))' * Vn(x+1:end-1,t+1) ...
                    + Vn(end,t+1)*(1.0-sum(pmf(1:(xmax-x),x+1,t+1)));
            end
            if ( x<xmax && t==tmax )
                %compute cost of no replacement
                notreplace = dot(pmf(1:(xmax-x),x+1,t),Vn(x+1:end-1,t)) ...
                + Vn(end,t)*(1.0-sum(pmf(1:(xmax-x),x+1,t)));
            end
            %use if so that you can compute thresholds in one go.
            if(x<xmax && replace<notreplace)
                Vn1(x+1,t) = replace;
                %set newly found threshold
                if (T(t)<1)
                    T(t)=x;
                end
            end
            if(x<xmax && replace>=notreplace)
                Vn1(x+1,t) = notreplace;
            end
        end
    end
    M = max(Vn1-Vn,[],"all");
    m = min(Vn1-Vn,[],"all");
    span = M-m;
    Vn = Vn1;
end

CostRate=(M+m)/2;
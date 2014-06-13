function [w,warm,gap]=lipm2(y,X,D,l1,l2,warm,verbose)
% function w=lipm2(y,X,D,l1,l2)
%
% A primal-dual interior point method for fused lasso
% y     target variable
% X     predictors
% D     sparse grouping matrix (one row per edge)
%           fused lasso  D_{k,i} = 1 D_{k,i+1} = -1
%           trend filter D_{k,i} = 2, D_{k,i-1} = -1, D_{k,i+1} = -1
% l1    penalty for singleton terms (vector) of size [p 1], p: row number
% in X
% l2    penalty for grouping/difference terms (vectors) of size [e 1], e:
% row number in D
%
% The problem being solved:
%    minimize    (1/2)||y - Xw||^2_2 + l1'||w||_1 + l2'||Dw||_1
%
% Reformulated problem:
%    minimize    (1/2)r'r + l1'*|z| + l2'*|d|
%    subject to  r = y - Xw
%                z = w
%                d = Dw
%
% Vladimir Jojic, 2009

N = size(X,1); assert(size(y,1) == N);
T = size(X,2);
E = size(D,1);

ALPHA     = 0.1;   % backtracking linesearch parameter (0,0.5]
BETA      = 0.9;    % backtracking linesearch parameter (0,1)
MU        = 2;      % IPM parameter: t update
MAXITER   = 100;     % IPM parameter: max iteration of IPM
MAXLSITER = 40;     % IPM parameter: max iteration of line search
TOL       = 1e-6;   % IPM parameter: tolerance

if nargin<7
    verbose = 0;
end
% m=2*size(X,1)+2*size(D,1);
m=size(X,1);

% dual variables for reformulated problem
rho = zeros(N,1);  % N residuals dim of target variable
mu = zeros(T,1);    % dual variable for (z=w)
nu = zeros(E,1);    % dual variable for (d=Dw)


% dual variables for dual of reformulated problem
eta1 = 1*ones(T,1);   % dual variable for (mu - l1) <= 0
eta2 = 1.001*ones(T,1);   % dual variable for (-mu - l1) <= 0
zeta1 = 1*ones(E,1);  % dual variable for (nu - Dl2) <= 0
zeta2 = 1.001*ones(E,1);  % dual variable for (-nu - Dl2) <= 0
xi = zeros(T,1); % dual variable for [X' I D']*[rho;mu;nu] = 0


if (~isempty(warm))
    rho = warm.rho;     nu = warm.nu;      mu = warm.mu;
    eta1 = warm.eta1;   eta2 = warm.eta2;
    zeta1 = warm.zeta1; zeta2 = warm.zeta2;
    xi = warm.xi;
end

t = 1e-10; pobj = Inf; dobj = 0;
step = Inf;
f1 = -mu - l1  ; f2 = mu -l1;
f3 = -nu - l2; f4 = nu -l2;

newmu = mu; newnu = nu;
newrho = rho;
neweta1 = eta1;   neweta2 = eta2;
newzeta1 = zeta1; newzeta2 = zeta2;
newxi = xi;

newf1 = f1;newf2 = f2;
newf3 = f3;newf4 = f4;



M = [speye(N,N)              sparse([],[],[],N,T,0)            sparse([],[],[],N,E,0)               X ;
    sparse([],[],[],T,N,0)  -(spdiag(eta1./f1 + eta2./f2))     sparse([],[],[],T,E,0)               -speye(T,T);
    sparse([],[],[],E,N,0)   sparse([],[],[],E,T,0)           -(spdiag(zeta1./f3 + zeta2./f4))      -D;
    X'                        -speye(T,T)                      -D'                                  sparse([],[],[],T,T,0);];


M0 = [sparse(N,N)              sparse([],[],[],N,T,0)            sparse([],[],[],N,E,0)               X ;
    sparse([],[],[],T,N,0)     sparse([],[],[],T,T,0)            sparse([],[],[],T,E,0)               -speye(T,T);
    sparse([],[],[],E,N,0)   sparse([],[],[],E,T,0)              sparse([],[],[],E,E,0)               -D;
    X'                        -speye(T,T)                      -D'                                    sparse([],[],[],T,T,0);];


gap = realmax;
if verbose,tic,end;
du = ones(size(M,1),1);


for iters=0:MAXITER
    if verbose,toc,end;
     dobj = -(1/2*rho'*rho - rho'*y);

%      pobj1 = (1/2*(y-X*xi)'*(y-X*xi)+(eta1+eta2)'*l1 + (zeta1+zeta2)'*l2);
%      pobj2 = (1/2*(y-X*xi)'*(y-X*xi)+(abs(xi))'*l1 + abs(D*xi)'*l2);
     pobj1 = 1/2*(y-X*xi)'*(y-X*xi);
     pobj2 = pobj1;
     pobj1 = pobj1 + (eta1+eta2)'*l1 + (zeta1+zeta2)'*l2;
     pobj2 = pobj2 + (abs(xi))'*l1 + abs(D*xi)'*l2;
     
     pobj = min(pobj1,pobj2);
     
     cgap = -f1'*eta1 - f2'*eta2 - f3'*zeta1 - f4'*zeta2;

     lastGap = gap;
     gap = min(pobj - dobj,cgap);

              
     if gap < TOL || (lastGap - gap) == 0
         if verbose
             fprintf('converged gap:%d\n',gap);
         end
         break;
     end

     if (step >= 0.2)
        t =max(2*m*MU/gap, 1.2*t);
     end

     dg = [ones(N,1); -(eta1./f1 + eta2./f2); -(zeta1./f3 + zeta2./f4); zeros(T,1)];
     M = M0 + spdiag(dg);
     
    r = -[+rho - y            + X*xi;
          1/t*(1./f1 - 1./f2) - xi;
          1/t*(1./f3 - 1./f4) - D*xi;
          X'*rho - mu - D'*nu];

    
      du = (M+1e-12*speye(size(M,1)))\r;
    
    drho = du(1:N);
    dmu =  du(N+1:N+T);
    dnu =  du(N+T+1:N+T+E);
    dxi =  du(N+T+E+1:end);

    deta1  = eta1./f1.*dmu - eta1 - (1/t)*(1./f1);
    deta2  = - eta2./f2.*dmu - eta2 - (1/t)*(1./f2);
    dzeta1 = zeta1./f3.*dnu - zeta1 - (1/t)*(1./f3);
    dzeta2 = - zeta2./f4.*dnu - zeta2 - (1/t)*(1./f4);


    rdual = -[ rho-y         + X*xi; 
              -eta1 + eta2   - xi; 
              -zeta1 + zeta2 - D*xi];
    rcent = -[-eta1.*(-mu - l1) - 1/t; 
              -eta2.*(mu - l1) - 1/t; 
              -zeta1.*(-nu - l2) - 1/t; 
              -zeta2.*(nu - l2) - 1/t;];
    rpri = -[X'*rho - mu - D'*nu];
    residual = [rdual;rcent;rpri];


    negIdx1 = (deta1<0);
    negIdx2 = (deta2<0);
    negIdx3 = (dzeta1<0);
    negIdx4 = (dzeta2<0);
    step = 1;

    if (any(negIdx1))
        step = min( step, 0.99*min(-eta1(negIdx1)./deta1(negIdx1)) );
    end

    if (any(negIdx2))
        step = min( step, 0.99*min(-eta2(negIdx2)./deta2(negIdx2)) );
    end

    if (any(negIdx3))
        step = min( step, 0.99*min(-zeta1(negIdx3)./dzeta1(negIdx3)) );
    end

    if (any(negIdx4))
        step = min( step, 0.99*min(-zeta2(negIdx4)./dzeta2(negIdx4)) );
    end

    for liter=1:MAXLSITER
        newrho = rho + step*drho;
         newmu = mu + step*dmu;
         newnu = nu + step*dnu;
         neweta1 = eta1 + step*deta1;       newzeta1 = zeta1 + step*dzeta1;
         neweta2 = eta2 + step*deta2;       newzeta2 = zeta2 + step*dzeta2;
         newxi = xi + step*dxi;

        newf1 = -newmu - l1;             newf3 = -newnu - l2;
        newf2 = newmu - l1;              newf4 = newnu - l2;


        rdual = -[newrho-y+X*newxi; 
                  -neweta1 + neweta2 - newxi; 
                  -newzeta1 + newzeta2 - D*newxi];
              
        rcent = -[-neweta1.*(-newmu - l1) - 1/t; 
                  -neweta2.*(newmu - l1) - 1/t; 
                  -newzeta1.*(-newnu - l2) - 1/t; 
                  -newzeta2.*(newnu - l2) - 1/t];
              
        rpri = -[X'*newrho - newmu - D'*newnu];
        newResidual = [rdual;rcent;rpri];


        if ( max(max(max(newf1),max(newf2)),max(max(newf3),max(newf4))) < 0 && ...
                norm(newResidual) <= (1-ALPHA*step)*norm(residual) )
            break;
        end

        step = BETA*step;
    end

    if verbose
        cost =     0.5*(y - X*newxi)'*(y-X*newxi) + l1'*abs(newxi) + l2'*abs(D*newxi);
        fprintf('it: %4d \t cost: %3.5e primal: %3.5e  dual:%3.5e gap: %3.5e  drho: %3.5e  dmu:  %3.5e  dnu: %3.5e \n',iters,cost,pobj,dobj,gap,sum(abs(newrho - rho)),sum(abs(mu-newmu)), sum(abs(nu-newnu)));
    end



    rho = newrho;    
     mu = newmu;  
     nu = newnu;
    eta1 = neweta1;      eta2 = neweta2;
    zeta1 = newzeta1;    zeta2 = newzeta2;
    xi = newxi;
    
    f1 = newf1;    f2 = newf2;
    f3 = newf3;    f4 = newf4;

end
w = xi;
warm.rho = rho;
warm.nu = nu;
warm.mu = mu;
warm.eta1 = eta1;
warm.eta2 = eta2;
warm.zeta1 = zeta1;
warm.zeta2 = zeta2;
warm.xi = xi;



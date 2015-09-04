%%(initializations)
% N = number of test cases
% error vector stored the estimation error for each of the N cases
% X is a set of N random points in the 2D plane
% description of variables already described in sdp.m has been omitted, as they are the same
function [error]=ml2(sigma,N)
error=zeros(1,N);
ml=zeros(1,N);
X=100*rand(N,2)-50;

% for-loop to estimate location for each of the test cases
for z = 1:1:N
    echo off
    clc
    clearvars -except z error X sigma 
    
    y=50*[1,-1, 1, -1 0 .25 .1 .2; ; 
       1, 1, -1, -1 0 .5 -.1 -.2];
    x1=X(z,:);           
    
    %taking the z-th point for estimation of location
    p0=-40; beta=4;d0=1;
    %sigma=4;
    rho=.5;
    
    %%
    %(initial calculations for L,lambda,d,P all calculated assuming known x1)
    % p0 is the power at unit length, taken to be 1m here
    % p = Vector of Recieved Signal Strenghts at each of the anchor nodes
    % d = Vector of distance of sensor in question from each of the anchor nodes
    % beta = path loss exponent, taken here to be 4dB
    % lambda is variable defined as 10^((P-P0)/10*beta).
    % q models the shadowing terms(shadowing is taken to be iid Gaussian
    % random variables). It is the covariance matrix of the same.
    % 
    p=p0*ones(1,8);
    d1=([x1' x1' x1' x1' x1' x1' x1' x1'] - y)'*([x1' x1' x1' x1' x1' x1' x1' x1'] - y); d1=sqrt((diag(d1))');
    p=p-10*beta*log10(d1/d0);
    lambda=10.^((p-p0*ones(1,8))/10*beta);
    L=diag(lambda);
    %%
    q=ones(8);
    for i=1:8
        for j=1:8
           if(i==j) q(i,j)=sigma^2;
           else q(i,j)=rho*sigma^2;
           end
        end
    end
    
    % The W is a weighting matrix which is proportional to the inverse of the covariance matrix
    W=((10*beta)^2 / (log(10))^2 )*inv(q);
    
    x=sym('x',[1,2]);
    for k=1:8
       tmp=x' ;
       tmp=tmp - y(:,k);  
       d(k)=(tmp(1)^2+tmp(2)^2)^0.5;  %d is convex, 4x1 vector of distances between x and y(i)
    end
    for j=1:8
       D(j,j)=[y(1,j),y(2,j),-1]*[eye(2) x';  x  x*x']*[y(1,j),y(2,j),-1]';
    end
    p=W*L*D*L';
    q=2*d'*ones(1,8)*W*L;
    F=trace(p)+trace(q);
    %% calculation of d vector, D matrix and final problem formulation in cvx(http://cvxr.com/cvx/doc/CVX.pdf) assuming x as a variable
    x=fminunc(@(x)F,0);
   
    error(z)=norm(x-x1); %storing the estimation error due to the z-th point
end

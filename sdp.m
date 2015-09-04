%%The following code implements the paper titled "RECEIVED SIGNAL STRENGTH-BASED SENSOR LOCALIZATION IN SPATIALLY CORRELATED SHADOWING" by vaghefi, buehrer. 
%Most studies for RSS (received signal strength) for sensor localization assume that the shadowing components are uncorrelated. However, here we assume 
%that the shadowing is spatially correlated. Under this condition, it can be shown that the localization accuracy can be improved if the correlation 
%among links is taken into consideration. The localization problem is formulated as an SDP minimization and solved using the cvx solver for
%MatLab. Just run to see results.

function [error]=sdp(sigma,N)
%%%%(initializations)                                                
%N=10;                                                  % N = number of test cases
error=zeros(1,N);                                       % localization error vector stored the estimation error for each of the N cases
X=100*rand(N,2)-50;                                     % X is a set of N random points (sensor locations to be estimated) in the 2D plane

%%
%%%% for-loop to estimate location for each of the test cases
for z = 1:1:N
    echo off                                            %allowing the commands to be viewed as they execute, for debugging
    %cvx_pause(false);
    %clc
    clearvars -except z error X  sigma  N          
    y=50*[1,-1,  1, -1, 0, .25,  .1,  .2;                      %8 anchors
          1, 1, -1, -1, 0,  .5, -.1, -.2];
    x1=X(z,:);                                          % randomly generated sensor point, to be estimated later
                                                        % taking the z-th point for estimation of location
    p0=-40;                                             % p0 is the power at unit length (d0=1)
    beta=4;                                             % beta = path loss exponent, taken here to be 4dB
    d0=1; 
    rho=.8;
    
    %%
    %%%%(initial calculations for L,lambda,d,P all calculated assuming known sensor x1)
    
    % Shadowing modeled as not-independant identically distributed Gaussian random
    % variables, with the covariance matrix q. For iid (independant
    % identically distributed Gaussian random variables, q=sigma^2 * eye(4)
  
    p=p0*ones(1,8);                                                             
    d1=(x1'*ones(1,8) - y)'*(x1'*ones(1,8) - y); d1=sqrt((diag(d1))');  % d1 =Vector of distance of sensor in question from each of the anchor nodes
    p=p-10*beta*log10(d1/d0);                                           % p = Vector of Recieved Signal Strenghts at each of the anchor nodes, (EQN #1)
    lambda=10.^((p-p0*ones(1,8))/10*beta);                              % lambda is a user-defined variable, defined as 10^((P-P0)/10*beta).
    L=diag(lambda);
 
    q=ones(8);                                                          % Filling the q matrix (EQN #2)
    for i=1:8
        for j=1:8
           if(i==j) q(i,j)=sigma^2;
           else q(i,j)=rho*sigma^2;
           end
        end
    end
    
    W=((10*beta)^2 / (log(10))^2 )*inv(q);                       % The W is a weighting matrix which is proportional to the inverse of the covariance matrix
    
    %% 
    %%%% calculation of d vector, D matrix and final problem formulation in cvx(http://cvxr.com/cvx/doc/CVX.pdf) assuming x as a variable
    
    %cvx_begin quiet sdp
    cvx_begin sdp
    variables x(1,2);
    variable pp(1,1);                                                       %pp  is z as in paper
    variable D(8,8);
    variable d(1,8);

                                                                       
    p=W*L*D*L';                                                             %p and q are cost functions to be minimized
    q=2*d'*ones(1,8)*W*L;
    minimize (trace(p) - trace(q))                                          %W- m x m; L- m x m; D - m x m; d - m x 1, m=no. of anchors

    subject to 
                                                                            %constraints
    [eye(2) x'; x pp] == semidefinite(3) ;                                  % 3rd constraint ( #14)
    
      for j=1:8                                                             %(eqn #13 in paper) 1st constraint
       D(j,j)>=[y(1,j),y(2,j),-1]*[eye(2) x';  x  pp]*[y(1,j),y(2,j),-1]';  %D matrix modeled to help making minimization problem linear and convex
      end 
      
    [D d'; d 1] == semidefinite(9) ;                                        %2nd constraint                               
    
    cvx_end
    
    error(z)=norm(x-x1);                                                    %storing the estimation error due to the z-th point
end

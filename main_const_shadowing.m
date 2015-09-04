                                                                            %main file with a fixed value of shadowing, sigma is a scalar
sigma=4; N=50;                                                             %takes less than a minute to run
error=sdp(sigma,N) ;                                                        %sdp estimator 
[a,b]=ecdf(error);                                                          %getting the CDF of the localization error
plot(b,'-ro')
xlabel('localization error')
ylabel('CDF')
hold                                                                        %ML estimator file, kept in same directory
error_ml=ml2(sigma,N);
[a,b]=ecdf(error_ml);
plot(b,'-.b');



%main file for varying shadowing, i.e. sigma as vector
%takes several minutes to run
sigma=1:5:41;
N=40;
rmse=zeros(1,length(sigma));                                                %rmse initialized to (1 x length(sigma)) to store root mean square values
for k=1:length(sigma)
error=sdp(sigma(k),N)                                                       %solve sdp for variable shadowing
rmse(k)=norm(error);                                                        %filling rmse values corresponding to values of shadowing
end
plot(sigma,rmse,'-ro')
xlabel('shadowing(DB)')
ylabel('RMSE')
hold
for k=1:length(sigma)
error=ml2(sigma(k),N)                                                       %solve for ml2 to get error vector and resultant rmse values 
rmse(k)=norm(error);
end
plot(sigma,rmse,'-.b')                                                      %plot on same graph

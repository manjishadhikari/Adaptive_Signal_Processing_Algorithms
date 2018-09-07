
function [R_normal, R_fb]=corr_matrix(X)

[M,L]=size(X);
R_normal=1/L*X*X';   %Normal Correlation Matrix

J=flipud(eye(M));
R_fb=1/(2*L)*(X*X'+J*conj(X)*transpose(X)*J);  %Correlation Matrix using Forward Backward Averaging 
return
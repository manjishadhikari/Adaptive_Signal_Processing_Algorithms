
% EECS 844 Exam-1
% Manjish Adhikari - 2870257
% Question 5

%%
clear;
close all;
load('V:\EECS-844\P5.mat');

M=15;   %Number of arrays
R=zeros(M,M);
for id=1:size(X,2)
  R=R+(X(:,id)*X(:,id)');
end
R=R/size(X,2);    % Correlation matrix R
eig_d0=eig(R);     %eigen values

max_eig_d0=max(abs(eig_d0));
min_eig_d0=min(abs(eig_d0));
cond_no_d0=max_eig_d0/min_eig_d0;


figure;plot(10*log10(eig_d0));title('Eigen values ');
xlabel('Number of eigen values');
ylabel('Eigen values in dB')


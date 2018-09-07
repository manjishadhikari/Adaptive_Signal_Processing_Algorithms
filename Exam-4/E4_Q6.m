%BIC and MUSIC using normal sample correlation matrix

clear all;
close all
load('V:\EECS-844\Exam-4\P6.mat');

for sample=1:4 
if sample==1
  X=X1;
elseif sample==2
  X=X2;
elseif sample==3
  X=X3;
else X=X4;
end
  
[M,N]=size(X);
p=M;      %Number of array elements
n=N;      %Number of time samples
R=1/n*(X*X');     %Correlation Matrix
eig_vals=eig(R);
eig_vals=sort(real(eig_vals),'descend'); %Sort if not sorted
% 
for q=1:p  
  term1=sum(log(eig_vals(q:p)));
  term2=(p-q+1).*log(sum(eig_vals(q:p)/(p-q+1)));
  Loglq=n*(term1-term2);
  BIC(q)=-2*Loglq+((q-1)*(2*p-q+1)+1)*log(n);
end

% figure;plot([1:p],(BIC)); title(sprintf('BIC for sample %d ',sample))
[min_val,min_idx]=min(BIC);

%Number of signals
p=min_idx-1;
fprintf('For sample %d no of signals is  %d\n',sample,p)

%MUSIC Implementation
[Q,D]=eig(R);
eig_values=real(diag(D));  %eigen values 
sorted_eig=(sort(eig_values,'descend')); %Sorted eigen values
if D(1)~=sorted_eig(1)
  Q=fliplr(Q);     %Arrange Q based on decreading eigen values 
end

sampling=600;  %Number of samples
phi=linspace(-pi/2,pi/2,sampling);
theta=pi*sin(phi);

for idx=1:sampling
  U=0;
  for k=p+1:M
    sv=transpose(exp((-1i)*theta(idx)*[0:M-1]));   %Steering vectors
    U=U+(abs(sv'*Q(:,k)))^2;
  end
  P(idx)=1/U;
end
figure(1);hold on; plot(theta *180/pi,20*log10(P));
end
figure(1); title('MUSIC Psedospectrum using Bayesian Information Criteria');
xlabel('Theta (deg)');
ylabel('Magnitude in dB')
legend('X1','X2','X3','X4')
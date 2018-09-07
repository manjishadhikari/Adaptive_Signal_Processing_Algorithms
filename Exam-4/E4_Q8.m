%BIC and MUSIC Implementation using Spatial Smoothing 

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
p=M;  %Number of antenna 
n=N;   %Number of time samples
Mbar=12; %Sub array size
K=M-Mbar+1; %Number of subarrays
Rss=zeros(Mbar,Mbar);
for subarrayidx=1:K
  
  Xnew=X(subarrayidx:subarrayidx+Mbar-1,:);  
  Rss=Rss+Xnew*Xnew';       %Correlation matrix using Spatial Smoothing
  
end
Rss=Rss/(n*K);

eig_vals=eig(Rss);
eig_vals=sort(real(eig_vals),'descend'); %Sort eig vals if not sorted
p=Mbar;
for q=1:p  
  term1=sum(log(eig_vals(q:p)));
  term2=(p-q+1)*log(sum(eig_vals(q:p)/(p-q+1)));
  Loglq=n*(term1-term2);
  BIC(q)=-2*Loglq+((q-1)*(2*p-q+1)+1)*log(n);
end

%figure;plot([1:p],(BIC))
[min_val,min_idx]=min(BIC);

%Number of signals
p=min_idx-1;
fprintf('For sample %d no of signals is  %d\n',sample,p)

%MUSIC Implementation
[Q,D]=eig(Rss);
eig_values=real(diag(D));  %eigen values 
sorted_eig=(sort(eig_values,'descend')); %Sorted eigen values
if D(1)~=sorted_eig(1)
  Q=fliplr(Q);     %Arrange Q based on decreasing eigen values 
end

sampling=600;
phi=linspace(-pi/2,pi/2,sampling);
theta=pi*sin(phi);

for idx=1:sampling
  U=0;
  for k=p+1:Mbar
    sv=transpose(exp((-1i)*theta(idx)*linspace(0,Mbar-1,Mbar)));   %Steering vectors
    U=U+(abs(sv'*Q(:,k)))^2;
  end
  P(idx)=1/U;
end
figure(1);hold on;plot(theta*180/pi,20*log10(abs(P)));
end
figure(1); title('MUSIC Psedospectrum using Spatial Smoothing');
xlabel('Theta (deg)');
ylabel('Magnitude in dB')
legend('X1','X2','X3','X4')

load('V:\EECS-844\Exam-4\P6.mat');
for idx=1:3
  if idx==1
    [p,R]=E4_Q6();
  elseif idx==2
    [p,R]=E4_Q7();
  else
    [p,R]=E4_Q8();
  end
  
%p=Number of signals

% fprintf('For sample %d no of signals is  %d\n',sample,p)
%MUSIC Implementation
[Q,D]=eig(R);
%Steering vectors
M=size(R,1)
sampling=600;
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
figure(100);hold on;plot(20*log10(P)); 
end
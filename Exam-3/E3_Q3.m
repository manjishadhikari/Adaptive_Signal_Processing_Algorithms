clear all
load('V:\EECS-844\Exam-3\P3.mat')
close all

N=10;      %Number of Antenna Elements
s1=zeros(N,1);    %Steering vector 1 at theta=-20
s2=zeros(N,1);     %Steering vector 2 at theta=+20
for k=1:N
  s1(k)=exp((-1i)*(k-1)*(-20*pi/180));
   s2(k)=exp((-1i)*(k-1)*(+20*pi/180));
end
w0=inv(R+s1*s1'+s2*s2')*s1;   %Optimum filter 

%% Beam Pattern Generation
phi=linspace(-pi/2,pi/2,14000);  %phi 
 theta=pi*sin(phi);       %theta

sv=zeros(N,length(phi));
for i=1:length(phi)
  for k=1:N
    sv(k,i)=exp(-1i*(k-1)*theta(i));  %Steering  vectors
  end
end
 
beam_pattern=sv'*w0;    %Beam Pattern 
figure(1);plot(theta*180/pi,20*log10(abs(beam_pattern)));
xlabel('THETA in degrees');ylabel('Gain in dB');
title('MVDR Beam Pattern')


%% Steepest Decent
eigen_vals=eig(R);
eig_max=max(eigen_vals);
mu_max=2/eig_max; 
mu=mu_max;   %Using mu_max
g=zeros(N,1);

num_iter=3000;           %Number of Iterations
w=zeros(10,num_iter+1);
for iter=1:num_iter
    g(:,iter)=2*(R*w(:,iter)+s1*s1'*w(:,iter)-s1+s2*s2'*w(:,iter));
    w(:,iter+1)=w(:,iter)-1/2*mu*g(:,iter);
    J(iter)=cost_func(R,s1,s2,w(:,iter+1));
    Dev(iter)=(w0-w(:,iter+1))'*(w0-w(:,iter+1));
end

figure(2);plot([1:length(Dev)],10*log10(Dev)); 
xlabel('Number of Iterations');ylabel('Squared dev(dB)')
title('Squared Deviation  in each iteration for mu\_max')
figure(5);plot([1:length(Dev)],10*log10(Dev));

%mu value greater than mu_max
mu=mu_max+0.00001;     
num_iter=3000;           %Number of Iterations
w=zeros(10,num_iter+1);
for iter=1:num_iter
    g(:,iter)=2*(R*w(:,iter)+s1*s1'*w(:,iter)-s1+s2*s2'*w(:,iter));
    w(:,iter+1)=w(:,iter)-1/2*mu*g(:,iter);
    J(iter)=cost_func(R,s1,s2,w(:,iter+1));
    Dev2(iter)=(w0-w(:,iter+1))'*(w0-w(:,iter+1));
end
figure(3);plot([1:length(Dev2)],10*log10(Dev2));
xlabel('Number of Iterations');ylabel('Squared dev(dB)')
title('Squared Deviation  in each iteration for mu>mu\_max')

%mu value less than mu_max
mu=mu_max-0.00001;      
num_iter=3000;           %Number of Iterations
w=zeros(10,num_iter+1);
for iter=1:num_iter
    g(:,iter)=2*(R*w(:,iter)+s1*s1'*w(:,iter)-s1+s2*s2'*w(:,iter));
    w(:,iter+1)=w(:,iter)-1/2*mu*g(:,iter);
    J(iter)=cost_func(R,s1,s2,w(:,iter+1));
    Dev3(iter)=(w0-w(:,iter+1))'*(w0-w(:,iter+1));
end
figure(4);hold on;plot([1:length(Dev3)],10*log10(Dev3));
xlabel('Number of Iterations');ylabel('Squared dev(dB)')
title('Squared Deviation  in each iteration for mu<mu\_max')

%% Steepest decent using optimum mu

num_iter=3000;           %Number of Iterations
w2=zeros(10,1);
for iter=1:num_iter
    g2=2*(R*w2+s1*s1'*w2-s1+s2*s2'*w2);
    mu_opt=(g2'*g2)/(g2'*(R+s1*s1'+s2*s2')*g2); %Optimum mu 

    w2=w2-1/2*mu_opt*g2;
    J2(iter)=cost_func(R,s1,s2,w2);
     Dev4(iter)=(w0-w2)'*(w0-w2);
end

figure(5);hold on;plot([1:length(Dev4)],10*log10(Dev4));
legend('w','w\_opt');title('Squared Deviation')
xlabel('Number of Iterations');ylabel('Squared Dev (dB)')

dev1_dB=10*log10(Dev);
dev4_dB=10*log10(Dev4);
idx1=find(dev1_dB<-30,1,'first');
idx2=find(dev4_dB<-30,1,'first');

disp(sprintf('Iterations For -30 dB squared deviation using mu_max: %d',idx1))
disp(sprintf('Iterations For -30 dB squared deviation using mu_opt: %d',idx2))
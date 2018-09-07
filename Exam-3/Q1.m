clear;
load('V:\EECS-844\Exam-3\P1.mat');

M=19;     %Length of each snapshot
K=length(d);
N=K-M+1;

X=complex(zeros(M,N));
D=complex(zeros(M,N));  
 for k=1:N
      X(:,k)=flipud(x(k:k+M-1));   %Snapshot matrix
      D(:,k)=flipud(d(k:k+M-1));   %Snapshot matrix
 end

 %% Using Linear Prediction
 
 R=1/N*D*D';     %Correlation Matrix
 
  r=zeros(M,1);    %Cross Correlation Matrix
    for i=M:K-1
      r=r+ flipud(d(i-M+1:i)).*conj(d(i+1));
    end
    r=r/N;
    
    wf=inv(R)*r;   %Weiner filter
    am=[1; -wf];   %Forward prediction error filter 
    
    
%% Inverse system ID using x as input and d as output
   R2=1/N*X*X';       %Correlation Matrix
   
   P=zeros(M,1);     %Cross Correlation Matrix
    for i=M:N
      P=P+ flipud(x(i-M+1:i)).*conj(d(i));
    end
    P=P/N;
   w=inv(R2)*P;     %Filter using Weiner Hopf Equation
   
   %% Plotting Filter Properties
   figure(1); 
   subplot(2,2,1); plot(abs(am)); title('Magnitude Response of LPE Filter')
   subplot(2,2,3);plot(angle(am)); title(' Phase of LPE Filter'); 
    subplot(2,2,2); plot(abs(w)); title('Magnitude Response of ISI Weiner Filter')
   subplot(2,2,4);plot(angle(w)); title(' Phase of ISI Weiner  Filter'); 
   figure(4); freqz(am,512,'whole'); title('Frequency Response of LPE Filter')
   figure(5); freqz(w,512,'whole'); title('Frequency Response of ISI Weiner Filter')
%Barttlett window
clear all;
close all;
load('V:\EECS-844\Exam-4\P4.mat');
K=[1 5 10 20];

for iter=1:4
  M=length(x)/K(iter);  %Numner of samples for PSD calc
  Xi=reshape(x,M,K(iter)); %Division into chunks
  
  for div=1:size(Xi,2)
    Xii=Xi(:,div);          %Samples for PSD calculation
    %  f=-1/2:0.001:1/2;
    f=linspace(0,1,size(x,1)/2);        %Frequency
    for idx=1:length(f)
      temp=0;
      for n=0:M-1
        temp=temp+Xii(n+1)*exp((-1i)*2*pi*f(idx)*n);
      end
      P(idx)=(abs(temp)^2)/length(Xii);  %Power SPectrum for each chunk
    end
    P_all(:,div)=P;     %Power Spectrum for all chunks
  end
  P_final=1/K(iter)*sum(P_all,2);  %Normalizing 
 
  w=2*pi*f;
  figure(1);hold on;plot(w,20*log10(abs(P_final)))
end
figure(1); title('Bartlett PSD ')
xlabel('Normalized Frequency rad/sample')
ylabel('Normalized Power in dB');
legend('K=1','K=5','K=10','K=20')
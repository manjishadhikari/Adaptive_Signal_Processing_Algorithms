

% EECS 844 Exam-1
% Manjish Adhikari - 2870257
% Question 6

%%
clear;
close all
load('V:\EECS-844\Exam-1\P6.mat');

M=14;   %Number of antenna array elements
 phi=linspace(-pi/2,pi/2,14000);  %phi 
 theta=pi*sin(phi);       %theta

sv=zeros(M,length(phi));
for i=1:length(phi)
  for k=1:M
    sv(k,i)=exp(-1i*(k-1)*theta(i));  %Steering  vector
  end
end
 
beam_pattern_non_adap=sv'*w_non_adap;    
beam_pattern_adap=sv'*w_adap;

figure(1);plot(theta*180/pi,20*log10(abs(beam_pattern_non_adap)),'r','DisplayName','Non Adaptive');
xlabel('THETA in degrees');ylabel('Gain in dB'); 
figure(2);plot(phi*180/pi,20*log10(abs(beam_pattern_non_adap)),'r','DisplayName','Non Adaptive');
xlabel('PHI in degrees');ylabel('Gain in dB'); 

figure(1);hold on;plot(theta*180/pi,20*log10(abs(beam_pattern_adap)),'b','DisplayName','Adaptive');title('theta vs gain ')
legend('show')
figure(2);hold on;plot(phi*180/pi,20*log10(abs(beam_pattern_adap)),'b','DisplayName','Adaptive');title('phi vs gain ')
legend('show')

figure(3);plot(theta*180/pi,(abs(beam_pattern_non_adap)),'r','DisplayName','Non Adaptive');
xlabel('THETA in degrees');ylabel('Gain ');
figure(4);plot(phi*180/pi,(abs(beam_pattern_non_adap)),'r','DisplayName','Non Adaptive');
xlabel('PHI in degrees');ylabel('Gain ');
figure(3);hold on;plot(theta*180/pi,(abs(beam_pattern_adap)),'b','DisplayName','Adaptive');title('theta vs gain ')
legend('show')
figure(4);hold on;plot(phi*180/pi,(abs(beam_pattern_adap)),'b','DisplayName','Adaptive');title('phi vs gain ')
legend('show')
load('V:\EECS-844\Exam-2\P3.mat')
close all;

L=length(X); %Length of snapshots
M=12;     %No of elements
N=30*M;   %Sampling 
sv=zeros(M,1);
power=zeros(N,1);
power_non_adap=zeros(N,1);

R=1/L*X*ctranspose(X); %Correlation matrix
R_inv=inv(R);

phi=linspace(-pi/2,pi/2,N);
theta=pi*sin(phi);

for i=1:N
  sv=transpose(exp(-j*theta(i)*[0:M-1]));   %Steering vectors
  for k=1:L
    power_non_adap(i)=power_non_adap(i)+(abs(ctranspose(sv)*X(:,k)))^2; %Power non adaptive
  end
    power_non_adap(i)=power_non_adap(i)/(L);
end
figure(1);plot(phi*180/pi,20*log10(abs(power_non_adap)));

sv_mvdr=zeros(N,1);
for i=1:N
  sv=transpose(exp(-j*theta(i)*[0:M-1]));
  sv_mvdr(i)=1./(ctranspose(sv)*R_inv*sv);  %MVDR power spectrum
end
figure(1);hold on;plot(phi*180/pi,20*log10(abs(sv_mvdr)));
legend('Non adaptive','MVDR')
xlabel('Spatial angle (phi in angles)');
ylabel('Power in dB')
title('MVDR and Non adaptive Power Spectrum')

phi_1=[-30 0 30]*pi/180;  %Three angles
angles=pi*sin(phi_1);    %Electrical Angle
w=zeros(M,3);
mvdr_beampattern=zeros(N,length(angles));
for i=1:length(angles)
  sv=transpose(exp(-j*angles(i)*[0:M-1]));
w(:,i)=R_inv*sv/(ctranspose(sv)*R_inv*sv);
%Beampattern
for l=1:N
  s=transpose(exp(-j*theta(l)*[0:M-1]));
  mvdr_beampattern(l,i)=ctranspose(s)*w(:,i);
end
end

figure(2);
for i=1:3
  subplot(3,1,i)
plot(phi*180/pi,20*log10(abs(mvdr_beampattern(:,i))));
title(sprintf('MVDR Beampattern for phi= %2d',round(phi_1(i)*180/pi)));xlabel('Theta in radians');ylabel('Power(dB)')
end
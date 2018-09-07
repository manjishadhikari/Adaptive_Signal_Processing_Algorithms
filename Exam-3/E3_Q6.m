load('V:\EECS-844\Exam-3\P6.mat')
close all;

% for X1
[M,~]=size(X1);
[R_normal1, R_fb1]=corr_matrix(X1);
eig_normal1=sort(eig(R_normal1));
eig_fb1=sort(eig(R_fb1));
[mvdr_beam_pattern_normal1]=Mvdr_power_spectrum(M,R_normal1);
[mvdr_beam_pattern_fb1]=Mvdr_power_spectrum(M,R_fb1);

phi=linspace(-pi/2,pi/2,600);
theta=pi*sin(phi);

figure(1);
subplot(5,1,1); plot(theta,10*log10(abs(mvdr_beam_pattern_normal1)));title('MVDR Power Spectrum for X1')
hold on; plot(theta,10*log10(abs(mvdr_beam_pattern_fb1)));xlabel('Theta');ylabel('Power (dB)')
legend('Normal','FB Averaging')

figure(2);subplot(5,1,1);plot([1:M],10*log10(abs(eig_normal1))); title('Eigen values for x1')
hold on; plot([1:M],10*log10(abs(eig_fb1))); ylabel('Eigen values in dB');xlabel('Eigen values Index')
legend('Normal','FB averaging')

% for X2
[M,~]=size(X2);
[R_normal2, R_fb2]=corr_matrix(X2);
eig_normal2=sort(eig(R_normal2));
eig_fb2=sort(eig(R_fb2));
[mvdr_beam_pattern_normal2]=Mvdr_power_spectrum(M,R_normal2);
[mvdr_beam_pattern_fb2]=Mvdr_power_spectrum(M,R_fb2);

figure(1);
subplot(5,1,2); plot(theta,10*log10(abs(mvdr_beam_pattern_normal2)));title('MVDR Power Spectrum for X2 ')
hold on; plot(theta,10*log10(abs(mvdr_beam_pattern_fb2)));xlabel('Theta');ylabel('Power (dB)')
legend('Normal','FB Averaging')

figure(2);subplot(5,1,2);plot([1:M],20*log10(abs(eig_normal2))); title('Eigen values for X2')
hold on;plot([1:M],20*log10(abs(eig_fb2))); ylabel('Eigen values in dB');xlabel('Eigen values Index')


% for X3
[M,~]=size(X3);
[R_normal3, R_fb3]=corr_matrix(X3);
eig_normal3=sort(eig(R_normal3));
eig_fb3=sort(eig(R_fb3));
[mvdr_beam_pattern_normal3]=Mvdr_power_spectrum(M,R_normal3);
[mvdr_beam_pattern_fb3]=Mvdr_power_spectrum(M,R_fb3);

figure(1);
subplot(5,1,3); plot(theta,10*log10(abs(mvdr_beam_pattern_normal3)));title('MVDR Power Spectrum for X3 ')
hold on; plot(theta,10*log10(abs(mvdr_beam_pattern_fb3)));xlabel('Theta');ylabel('Power (dB)')
legend('Normal','FB Averaging')

figure(2);subplot(5,1,3);plot([1:M],20*log10(abs(eig_normal3))); title('Eigen values for X3')
hold on;plot([1:M],20*log10(abs(eig_fb3)));  ylabel('Eigen values in dB');xlabel('Eigen values Index')


% for X4
[M,~]=size(X4);
[R_normal4, R_fb4]=corr_matrix(X4);
eig_normal4=sort(eig(R_normal4));
eig_fb4=sort(eig(R_fb4));
[mvdr_beam_pattern_normal4]=Mvdr_power_spectrum(M,R_normal4);
[mvdr_beam_pattern_fb4]=Mvdr_power_spectrum(M,R_fb4);

figure(1);
subplot(5,1,4); plot(theta,10*log10(abs(mvdr_beam_pattern_normal4)));title('MVDR Power Spectrum for X4')
hold on; plot(theta,10*log10(abs(mvdr_beam_pattern_fb4)));xlabel('Theta');ylabel('Power (dB)')
legend('Normal','FB Averaging')

figure(2);subplot(5,1,4);plot([1:M],20*log10(abs(eig_normal4))); title('Eigen values for X4')
hold on;plot([1:M],20*log10(abs(eig_fb4)));  ylabel('Eigen values in dB');xlabel('Eigen values Index')

% for X5
[M,~]=size(X5);
[R_normal5, R_fb5]=corr_matrix(X5);
eig_normal5=sort(eig(R_normal5));
eig_fb5=sort(eig(R_fb5));
[mvdr_beam_pattern_normal5]=Mvdr_power_spectrum(M,R_normal5);
[mvdr_beam_pattern_fb5]=Mvdr_power_spectrum(M,R_fb5);

figure(1);
subplot(5,1,5); plot(theta,10*log10(abs(mvdr_beam_pattern_normal5)));title('MVDR Power Spectrum for X5')
hold on; plot(theta,10*log10(abs(mvdr_beam_pattern_fb5)));xlabel('Theta');ylabel('Power (dB)')
legend('Normal','FB Averaging')

figure(2);subplot(5,1,5);plot([1:M],20*log10(abs(eig_normal5))); title('Eigen values for X5')
hold on;plot([1:M],20*log10(abs(eig_fb5)));  ylabel('Eigen values in dB');xlabel('Eigen values Index')





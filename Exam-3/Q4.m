clear all
 load('V:\EECS-844\Exam-3\P4.mat')
M = length(x);
mf = conj(flipud(x))./(x'*x);
K = 3;
mmfLen = M*K;
mmfOff = 2*M;
A = toeplitz([x; zeros(mmfLen-1,1)],[x(1) zeros(1,mmfLen-1)]);
Az = A;
Az(mmfOff-1,:) = 0;
Az(mmfOff+1,:) = 0;
[NA MA] = size(A);
[V D] = eig(A'*A);
lam_max = max(diag(D));
Phi = inv(A'*A)*A';
Phiz = inv(Az'*Az)*Az';
Phid = inv(A'*A + (0.02.*lam_max).*eye(MA))*A';
Phidz = inv(Az'*Az + (0.02.*lam_max).*eye(MA))*Az';
e = zeros(NA,1);
e(mmfOff) = 1;
mmf = Phi*e;
mmf = mmf./sqrt(mmf'*mmf)./sqrt(x'*x);
mmfz = Phiz*e;
mmfz = mmfz./sqrt(mmfz'*mmfz)./sqrt(x'*x);
mmfd = Phid*e;
mmfd = mmfd./sqrt(mmfd'*mmfd)./sqrt(x'*x);
mmfdz = Phidz*e;
mmfdz = mmfdz./sqrt(mmfdz'*mmfdz)./sqrt(x'*x);
mf_out = conv(mf,x);
mmf_out = conv(mmf,x);
mmfz_out = conv(mmfz,x);
mmfd_out = conv(mmfd,x);
mmfdz_out = conv(mmfdz,x);
MMloss = -20*log10(max(abs(mmf_out))/max(abs(mf_out)))
MMlossz = -20*log10(max(abs(mmfz_out))/max(abs(mf_out)))
MMlossd = -20*log10(max(abs(mmfd_out))/max(abs(mf_out)))
MMlossdz = -20*log10(max(abs(mmfdz_out))/max(abs(mf_out)))
yf_mf = conv(mf,y);
ymf_shft = [zeros(1*M-1,1); yf_mf];
yf_mmf = conv(mmf,y);
yf_mmfz = conv(mmfz,y);
yf_mmfd = conv(mmfd,y);
yf_mmfdz = conv(mmfdz,y);
len = length(y);
figure(1)
subplot(5,1,1)
plot(-(M-1):(M-1),20*log10(abs(mf_out)))
axis([-((K+1)*M/2-1) ((K+1)*M/2-1) -50 5])
grid on
title('matched filter')
xlabel('delay')
ylabel('filter response (dB)')
subplot(5,1,2)
plot(-((K+1)*M/2-1):((K+1)*M/2-1),20*log10(abs(mmf_out)))
axis([-((K+1)*M/2-1) ((K+1)*M/2-1) -50 5])
title('LS mismatched filter')
grid on
xlabel('delay')
ylabel('filter response (dB)')
subplot(5,1,3)
plot(-((K+1)*M/2-1):((K+1)*M/2-1),20*log10(abs(mmfz_out)))
axis([-((K+1)*M/2-1) ((K+1)*M/2-1) -50 5])
title('LS mismatched filter with zeroed rows')
grid on
xlabel('delay')
ylabel('filter response (dB)')
subplot(5,1,4)
plot(-((K+1)*M/2-1):((K+1)*M/2-1),20*log10(abs(mmfd_out)))
axis([-((K+1)*M/2-1) ((K+1)*M/2-1) -50 5])
title('LS mismatched filter with diagonal loading')
grid on
xlabel('delay')
ylabel('filter response (dB)')
subplot(5,1,5)
plot(-((K+1)*M/2-1):((K+1)*M/2-1),20*log10(abs(mmfdz_out)))
axis([-((K+1)*M/2-1) ((K+1)*M/2-1) -50 5])
title('LS mismatched filter with zeroed rows & diagonal loading')
grid on
xlabel('delay')
ylabel('filter response (dB)')
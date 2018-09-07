function [mvdr_beam_pattern]=Mvdr_power_spectrum(M,R)

phi=linspace(-pi/2,pi/2,600);
theta=pi*sin(phi);
for j=1:length(theta)
 for k=1:M
  sv(k,j)=transpose(exp((-1i)*theta(j)*(k-1)));
 end
  mvdr_beam_pattern(j)=1/((sv(:,j)'*inv(R)*sv(:,j)));
  
end
 mvdr_beam_pattern_dB=10*log10( mvdr_beam_pattern);
end

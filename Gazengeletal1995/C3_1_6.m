%
% 3.1.6 Examples

clear
%% Pipe parameter
L = 0.343/2;
c = 343;
tau = L/c;

%%



%% Sampling parameters
fs = 11e3;
fs = 10.5e3;
% Fmax = 5e3;
% fs = Fmax*2;
Fmax = fs/2;
dt = 1/fs;
NFFT = 10;
N = 2^NFFT; 
f = (0:N/2-1)/(N/2) * fs/2;
omega = 2*pi*f;

%==== Reflection coefficients
Rp = -exp(-2j*omega*tau); %% 

%==== Windowing / Low pass
[~, i] = min(abs(f-Fmax));
w = hann(i*2).';
Rp(1:i) = Rp(1:i) .* w(end/2+1:end);
Rp(i+1:end) = 0;


%==== complete it with its conjugate symmetry
Rp = [Rp flip(conj(Rp))]; 

%==== Plot 
figure(1)
ax(1) = subplot(2,1,1);
plot([f f+fs/2], real(Rp));
ylabel('\Re(R_p)', 'fontsize', 20)
ax(2) = subplot(2,1,2);
plot([f f+fs/2], imag(Rp))
ylabel('\Im(R_p)', 'fontsize', 20)
xlabel('Frequency (Hz)')
linkaxes(ax, 'x')
xlim([0, fs])


%% IFFT
figure(2)
rd = ifft(Rp,  'symmetric');

%==== low pass filter
% Rp = 
% rd = lowpass(rd, 2000, fs,'ImpulseResponse','iir','Steepness',0.5);
sum(rd(1:end/2))

%==== Plot
t = (0:length(rd)-1)*1/fs;
plot(t, rd, 'LineWidth', 2)
hold on;
%%
t = 0:1e-6:0.01;
rp = zeros(length(t), 1);
rp(ceil(2*tau/1e-6)) = -1;

plot(t, rp, '--', 'LineWidth', 2)

xlim([0 0.01])


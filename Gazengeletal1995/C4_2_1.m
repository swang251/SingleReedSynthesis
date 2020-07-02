%Implementation of the example in Gazengel et al. (1995) Time Domain
%Simulation of Single Reed Instrument.
% C4.2.1 Numerical Simulation Using a Lossless Cylindrical Tube
%
% by Song Wang, McGill University, Jun. 2020

clear
%% Parameters
%==== Physical Parameters
rho = 1.18;
c = 343;
H = 0.4e-3;
r = 0.0075;
w = 0.0125;
S = pi*r^2;
Zc = rho*c/S;

%==== Reed parameters
mur = 0.0231;   % kg/m2
omegar = 20000; % rad/s
kr = mur*omegar^2;
Pm = 1500;      % pa
Pc = Pm - kr*H;

%==== Pipe parameters
L = 0.343/2;
c = 343;
tau = L/c;

%% Reflection function calculation
%==== Sampling Parameters
fs = 10.5e3; Fmax = fs/2;
fs = 11e3; Fmax = fs/2;
% fs = 48e3; Fmax = fs/2;
% fs = 48e3; Fmax = 5e3;
dt = 1/fs;
NFFT = 15;
N = 2^NFFT; 
f = (0:N/2-1)/(N/2) * fs/2;
omega = 2*pi*f;
%==== Reflection coefficients
Rp = -exp(-2j*omega*tau); %% Eq. 19

%==== Windowing / Low pass
[~, i] = min(abs(f-Fmax));
Rp(i+1:end) = 0;  % Zero-padding
% win = hann(i*2).';
% Rp(1:i) = Rp(1:i) .* win(end/2+1:end);    % windowing

%==== complete it with its conjugate symmetry
Rp = [Rp flip(conj(Rp))];
    %==== Plot 
    figure(1)
    ax(1) = subplot(5,1,1);
    plot([f f+fs/2], real(Rp));
    xline(fs/2, 'r', 'LineWidth', 2);
    ylabel('\Re(R_p)', 'fontsize', 20)
    xlabel('Frequency (Hz)', 'fontsize', 20)
    
    ax(2) = subplot(5,1,2);
    plot([f f+fs/2], imag(Rp))
    xline(fs/2, 'r', 'LineWidth', 2);
    ylabel('\Im(R_p)', 'fontsize', 20)
    xlabel('Frequency (Hz)', 'fontsize', 20)
    linkaxes(ax, 'x')
    xlim([0, fs])

%==== IFFT
r = ifft(Rp,  'symmetric');
% r = lowpass(r, Fmax, fs,'ImpulseResponse','iir','Steepness',0.5);
t = (0:length(r)-1)*1/fs;
[~, trun] = max(t(t<10*tau));
sum(r(1:trun))
r = r(1:trun);
t = t(1:trun);
    %==== Plot
    subplot(5,1,3);
    plot(t, r, 'LineWidth', 2)
    ylabel('r', 'fontsize', 20)
    xlabel('t(s)', 'fontsize', 20)
    hold on;


%% Simulation
time = 0.1;
samples = floor(time * fs);
p = zeros(1, samples);
U = zeros(1, samples);
y = zeros(1, samples);
ppast = 0;
pn_ = 0;
NewtonRaphson = false;
%==== Main Loop
tic
for n = 2:samples
    iconv = min([length(r), n-1]); 
    ppast = sum( r(1:iconv) .* (p(n-1:-1:n-iconv)+Zc*U(n-1:-1:n-iconv)) ); % this is p_past[n+1]
    
    p(n) = ppast + Zc*U(n-1);
    if abs(p(n)/Pm)<0.7 && NewtonRaphson
        flag = true;
        while (abs(pn_ - p(n))/pn_ > 1e-3) || flag
            flag = false;
            sgn = sign(Pm-p(n));
            dp = ( sgn*0.5*1/sqrt(abs(Pm-p(n)))*(p(n)-Pc) + sqrt(abs(Pm-p(n))) )...
                      * Zc * sqrt(2/rho)*w / kr * sgn - 1;
            if abs(dp) > 20
                break;
            end
            y(n) = (p(n)-Pm)/kr;  % Eq. (53)
            if y(n) <= -H % closed reed
                U(n) = 0;
            else
                U(n) = sqrt(2/rho) * w * sqrt(abs(Pm-p(n))) * (y(n)+H) * sign(Pm-p(n));
        %         U(n) = sqrt(2/rho) * w/mur/omegar^2 * sqrt(abs(Pm-p(n))) * (p(n)-Pc) * sign(Pm-p(n)); % Eq. 54
            end
            pn_ = p(n);
            p(n) = pn_ - (ppast + Zc*U(n) - pn_) / dp;
        end
    end
    y(n) = (p(n)-Pm)/kr;  % Eq. (53)
    if y(n) <= -H % closed reed
        U(n) = 0;
    else
        U(n) = sqrt(2/rho) * w * sqrt(abs(Pm-p(n))) * (y(n)+H) * sign(Pm-p(n));
%         U(n) = sqrt(2/rho) * w/kr * sqrt(abs(Pm-p(n))) * (p(n)-Pc) * sign(Pm-p(n)); % Eq. 54
    end    
end
toc
%==== plot
figure(1)
subplot(5,1,4)
title(['Pm = ' num2str(Pm)])
plot((0:dt:time-dt), p/Pm); hold on;
ylim([-2,2])
ylabel('p/P_m', 'fontsize', 20)
xlabel('t(s)', 'fontsize', 20)
subplot(5,1,5)
plot((0:dt:time-dt), y/H); hold on;
ylabel('y/H', 'fontsize', 20)
xlabel('t(s)', 'fontsize', 20)

% sound(p,fs)
%% Newton-Raphson Method


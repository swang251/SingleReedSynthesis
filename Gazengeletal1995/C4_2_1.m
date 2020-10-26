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

type = input("1) incontinuous, 2) continuous, 3) low-passed: ");
%% Reflection function calculation
%==== Sampling Parameters
if type == 1 || type == 3
    fs = 10.5e3; 
    disp("Fs = 10.5 kHz")
elseif type == 2
    fs = 11e3; 
    disp("Fs = 11 kHz")
else
    error("please enter the right type!");
end
% fs = 44000; 
% fs = 48e3; 
% fs = 96e3; 
% fs = 48e3; Fmax = 5e3;
Fmax = fs/2;
dt = 1/fs;
NFFT = 15;
N = 2^NFFT; 
f = (0:N/2-1)/(N/2) * fs/2;
omega = 2*pi*f;
%==== Reflection coefficients
Rp = -exp(-2j*omega*tau); %% Eq. 19

if type == 3
%==== Windowing / Low pass
    [~, i] = min(abs(f-Fmax));
    Rp(i+1:end) = 0;  % Zero-padding
    win = hann(i*2).';
    Rp(1:i) = Rp(1:i) .* win(end/2+1:end);    % windowing
end

%==== complete it with its conjugate symmetry
Rp = [Rp flip(conj(Rp))];
    %==== Plot 
    figure(1)
    ax(1) = subplot(3,1,1);
    plot([f f+fs/2], real(Rp));
    xline(fs/2, 'r', 'LineWidth', 2);
    ylabel('\Re(R_p)', 'fontsize', 20)
    xlabel('Frequency (Hz)', 'fontsize', 20)
    
    ax(2) = subplot(3,1,2);
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
    subplot(3,1,3);
    plot(t, r, 'LineWidth', 2)
    ylabel('r', 'fontsize', 20)
    xlabel('t (s)', 'fontsize', 20)
    hold on;

%% Simulation
time = 1;
samples = floor(time * fs);
p = zeros(1, samples);
U = zeros(1, samples);
y = zeros(1, samples);
ppast = 0;
pn_ = 0;
%==== Main Loop

tic
for n = 2:samples
    iconv = min([length(r), n-1]); 
    ppast = sum( r(1:iconv) .* (p(n-1:-1:n-iconv)+Zc*U(n-1:-1:n-iconv)) ); % this is p_past[n+1]
        
    pn = p(n-1);
%     flag = true;
%     while (abs(pn_ - pn)/pn_ > 1e-5) || flag
%         flag = false;
    %===== Newton-Raphson's method
    for iter = 1:50
        y(n) = (pn-Pm)/kr;  % Eq. (53)
        if y(n) <= -H % closed reed
            U(n) = 0;
        else
            U(n) = sqrt(2/rho) * w * sqrt(abs(Pm-pn)) * (y(n)+H) * sign(Pm-pn);
%                 U(n) = sqrt(2/rho) * w/mur/omegar^2 * sqrt(abs(Pm-pn)) * (pn-Pc) * sign(Pm-pn); % Eq. 54
        end
        sgn = sign(Pm-pn);
        fp_deri = (-sgn*0.5/sqrt(abs(Pm-pn))*(pn-Pc) + sqrt(abs(Pm-pn)) )...
                  * Zc * sqrt(2/rho)*w / kr * sgn - 1;
        if abs(fp_deri) > 1000000 % pn ~ Pm
%         if abs(Pm-pn) / Pm < 0.00001
            break;
        end
        pn_ = pn;
        fp = (ppast + Zc*U(n) - pn_);
        pn = pn_ -  fp / fp_deri;
    end
       
%     pn = fzero(@(x) ppast + Zc * sqrt(2/rho)*w/kr*sqrt(abs(Pm-x))*(x - Pc)*sign(Pm - x)-x, p(n-1) );
    
    p(n) = pn;    
    y(n) = (pn-Pm)/kr;  % Eq. (53)
    if y(n) <= -H % closed reed
        U(n) = 0;
    else
        U(n) = sqrt(2/rho) * w/mur/omegar^2 * sqrt(abs(Pm-pn)) * (pn-Pc) * sign(Pm-pn); % Eq. 54
    end
    
  
end
toc
%==== plot
figure(2)
subplot(2,1,1)
title(['Pm = ' num2str(Pm)])
plot((0:dt:time-dt), p/Pm); hold on;
% ylim([-2,2])
ylabel('p/P_m', 'fontsize', 20)
xlabel('t (s)', 'fontsize', 20)
xlim([0, 0.04])
subplot(2,1,2)
plot((0:dt:time-dt), y/H); hold on;
ylabel('y/H', 'fontsize', 20)
xlabel('t (s)', 'fontsize', 20)
xlim([0, 0.04])

%%
sound(p,fs)
disp("synthesized fundamental frequency: " + num2str(mean(pitch(p.', fs, 'Range', [50, 1000]))))
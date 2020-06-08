%Clarinet Synthesis
% Implementing the method proposed by McIntyre, M.E. and Schumacher, R.T. 
% (1983) On the oscillation of musical instrument
%
% by Song Wang, McGill University, Jun. 2020
clear

L = 0.5;            % Pipe Length
c = 343;            % Sound speed
rho = 1.18;         % air density
T = 2 * L / c;      % round-trip time
Zc = 1;             % characteristic impedance (normalized)

fs = floor(128/T);
ts = 1/fs;

percent = 0.2;      
b = (percent * T / log(2)^(1/2) / 2) ^(-2); % Derived from Eq. (19)

t = 0:ts:T*2;     % length of impulse response equals the period
rEXP = exp(-b * (t-T).^2);

a = -1 / trapz(rEXP); % integration using the trapezoidal rule
r = a * rEXP;       % Eq. (18)
figure

plot(t/T, r);
xlabel('t/T');
ylabel('r');
%%
k = 0.2;            % parameters in Eq. (20)
qc = -2;
p = 3;
time = 0.05;
samples = floor(time * fs);
%==== Initialization
f = zeros(1, samples);
q = zeros(1, samples);
qo = zeros(1, samples);
qi = zeros(1, samples);
qh = 0;
%==== Main Loop
for i = 2:time*fs
    iconv = min([length(r), i-1]);
    qh = sum( r(1:iconv) .* (q(i-1:-1:i-iconv) + Zc*f(i-1:-1:i-iconv)));   % Eq. (11)
    q(i) = qh + Zc*f(i-1);    % Eq. (10)
    qo(i) = (q(i) + Zc * f(i))/2;
    qi(i) = (q(i) - Zc * f(i))/2;
    if q(i)<qc      % Eq. (1) + (20) with the discontinuity
        f(i)=0;
    else
        f(i) = k*(p-q(i))*(q(i)-qc);
    end
end

%% Plot
figure
subplot(4,1,1);
plot((1:samples)/fs, q)
subplot(4,1,2);
plot((1:samples)/fs, f)
subplot(4,1,3);
plot((1:samples)/fs, qo)
subplot(4,1,4);
plot((1:samples)/fs, qi)

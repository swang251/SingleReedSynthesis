%Clarinet Synthesis
% Implementing the method proposed by McIntyre, M.E., Schumacher, R.T. and
% Woodhouse, J. (1983) On the oscillation of musical instruments
%
% Variables:
%   f: flow
%   q: mouthpiece pressure
%   qi: incoming pressure
%   qo: outgoing pressure
%   Z: input impedance
%   p: mouth pressure
%   qc: extinction threshold
%
% The idealized clarinet is modeled with four equations
%   f   = F(q)                      Eq. (1)
%   qi  = r ** qo                   Eq. (4), ** means convolution
%   q   = qo + qi                   Eq. (7)
%   Z*f = qo - qi                   Eq. (8)
%
% Written in a compact form
%   f   = F(q)                      Eq. (1)
%   q   = qh + Z*f                  Eq. (10)
%   qh  = r ** {q + Z*f}            Eq. (11)
% 
% II. B Clarinet-like oscillations example
%   F(q) = k (p - q) ( q - qc),     for q>=qc,      Eq. (20)
%          0                  ,     for q <qc
% Combinining Eq. (10) and Eq. (20)
%   q(t) = qh + Z*k*(p-q)*(q-qc),   for qh >= qc, 
%          qh                   ,   for qh <  qc.
% Rewrriten in quadratic form
%   q^2 + Bq + C = 0,
%   where B = qc - p + 1/(Z*k), C = p*qc - qh/(Z*k)
% Solve the equation explicitly
%   q = (-B + sqrt(B^2 - 4C)) / 2,  for qh >= qc, or B^2-4C >=0
%       qh                       ,  for qh <  qc, or B^2-4C < 0
%
% Calculation Steps:
%   1. Update qh[n] using Eq. (11) with history values, [n-1, n-2,...]
%   2. Solve q[n] explicitly
%   3. Calculate f[n], qi[n] and qo[n]
%
% by Song Wang, McGill University, Jun. 2020
clear

L = 0.5;            % Pipe Length
c = 343;            % Sound speed
rho = 1.18;         % air density
T = 2 * L / c;      % round-trip time
Z = 1;             % characteristic impedance (normalized)

fs = floor(128/T);
ts = 1/fs;

percent = 0.05;     % Fig. 5(a)
percent = 0.2;      % Fig. 5(b)
% percent = 0.005;      % extra testing
b = (percent * T / log(2)^(1/2) / 2) ^(-2); % Derived from Eq. (19)

t = 0:ts:T*2;     % length of impulse response equals the period
rEXP = exp(-b * (t-T).^2);

a = -1 / trapz(rEXP);   % integration using the trapezoidal rule
r = a * rEXP;           % Eq. (18)

figure(1)
plot(t/T, r);
xlabel('t/T');
ylabel('r');
%%
k = 0.2;            % parameters in Eq. (20)
qc = -2;
p = 3;
B = -qc - p + 1/(Z*k);
time = 0.05;
samples = floor(time * fs);

%==== Initialization
f = zeros(1, samples);
q = zeros(1, samples);
qo = zeros(1, samples);
qi = zeros(1, samples);
qh = 0;

%==== Main Loop
tic
for i = 2:time*fs
    iconv = min([length(r), i-1]);
    qh = sum( r(1:iconv) .* (q(i-1:-1:i-iconv) + Z*f(i-1:-1:i-iconv)));   % Eq. (11)
    C = p*qc - qh/(Z*k);
    DELTA = B^2 - 4*C;
    if DELTA < 0
        q(i) = qh;
    else
        q(i) = (-B + sqrt(DELTA)) / 2;  % The solution of Eq. (1) + (20) with the discontinuity
    end
    qo(i) = (q(i) + Z * f(i))/2;
    qi(i) = (q(i) - Z * f(i))/2;
    if q(i)<qc      % Eq. (1) + (20) with the discontinuity
        f(i)=0;
    else
        f(i) = k*(p-q(i))*(q(i)-qc);
    end
    qi(i) = q(i) - qo(i);
end
toc

%==== Plot
figure(2)
subplot(4,1,1);
plot((1:samples)/fs, q); hold on;
ylabel('q')
subplot(4,1,2);
plot((1:samples)/fs, f); hold on;
ylabel('f')
subplot(4,1,3);
plot((1:samples)/fs, qo); hold on;
ylabel('q_o')
subplot(4,1,4);
plot((1:samples)/fs, qi); hold on;
ylabel('q_i')
xlabel('t(s)')

set(findall(gcf,'-property','FontSize'),'FontSize',18)
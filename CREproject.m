clc; clear;

% Constants
R = 8.314;       % J/mol-K
rho = 1000;      % mol/m³
Cp = 120;        % J/mol-K
UA = 0;
Tcool = 600;

% Heats of reactions (J/mol)
dH1 = -5000; dH2 = 120000; dH3 = 104000; dH4 = 30000;

% Rate constants (1/s)
k1 = 0.5; k2 = 0.2; k3 = 0.1; k4 = 0.3;

% Catalyst deactivation
kd = 0.01; % 1/m

% Initial conditions
Pn0 = 1.0; Pi0 = 0.0; N0 = 0.5; A0 = 0.0; G0 = 0.0; T0 = 700; a0 = 1.0;
y0 = [Pn0; Pi0; N0; A0; G0; T0; a0];
zspan = [0 10];

% ODE system
reform_model = @(z, y) [
    -k1*y(7)*y(1) - k2*y(7)*y(1) - k4*y(7)*y(1);           % dPn
     k1*y(7)*y(1);                                         % dPi
    -k3*y(7)*y(3);                                         % dN
     k2*y(7)*y(1) + k3*y(7)*y(3);                          % dA
     k4*y(7)*y(1);                                         % dG
    (-(dH1*k1*y(7)*y(1) + dH2*k2*y(7)*y(1) + ...
       dH3*k3*y(7)*y(3) + dH4*k4*y(7)*y(1))) / (rho*Cp);   % dT
    -kd * y(7)                                             % da
];

% Solve
[z, y] = ode45(reform_model, zspan, y0);
Pn = y(:,1); Pi = y(:,2); N = y(:,3); A = y(:,4); G = y(:,5); T = y(:,6); a = y(:,7);

% Derived: Conversion and Yield
Xpn = (Pn0 - Pn) / Pn0;
Ya = A / Pn0;

% Figure 1: Mole flows
figure(1)
plot(z, [Pn Pi N A G], 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Mole Flow (mol/s)')
legend('Pn','Pi','N','A','G')
title('Mole Flow Profiles'), grid on

% Figure 2: Temperature
figure(2)
plot(z, T, 'r', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Temperature (K)')
title('Temperature Profile'), grid on

% Figure 3: Catalyst activity
figure(3)
plot(z, a, 'k', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Catalyst Activity')
title('Deactivation Profile'), grid on

% Figure 4: Conversion and Yield
figure(4)
plot(z, Xpn, 'b--', z, Ya, 'm-', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Fraction')
legend('Conversion of Pn', 'Yield of Aromatics')
title('Conversion & Yield'), grid on

% Figure 5: All in one
figure(5)

subplot(2,2,1)
plot(z, [Pn Pi N A G], 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Mole Flow (mol/s)')
legend('Pn','Pi','N','A','G')
title('Mole Flow Profiles'), grid on

subplot(2,2,2)
plot(z, T, 'r', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Temperature (K)')
title('Temperature Profile'), grid on

subplot(2,2,3)
plot(z, a, 'k', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Catalyst Activity')
title('Deactivation Profile'), grid on

subplot(2,2,4)
plot(z, Xpn, 'b--', z, Ya, 'm-', 'LineWidth', 2)
xlabel('Reactor Length (m)'), ylabel('Fraction')
legend('Conversion of Pn', 'Yield of Aromatics')
title('Conversion & Yield'), grid on

sgtitle('Catalytic Reforming of Naphtha – All Graphs')

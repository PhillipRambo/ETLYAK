% Loudspeaker Cabinet Simulator
% Clearing the workspace
clear all;

% Fundamental physical constants
rho = 1.18; % Air mass density (kg/m^3)
c = 345; % Speed of sound (m/s)
pREF = 20e-6; % Reference sound pressure (Pa)

% Input parameters
UG = 2.83; % Voltage across terminals (V)
f = 10:1000; % Frequency scale (Hz)
s = 1i*2*pi*f; % Laplace operator (rad/s) using '1i' for the imaginary unit

% Loudspeaker parameters
RE = 5.6; % Coil DC resistance (Ohm)
LE = 0.65e-3; % Coil self-inductance (H)
BL = 5.1; % Force factor (N/A)
MMS = 0.007; % Diaphragm mass (kg)
CMS = 1.45e-3; % Suspension compliance (m/N)
RMS = 0.6; % Mechanical losses (Ns/m)
SD = 0.0054; % Diaphragm area (m^2)

% Simulation starts here
% Use a switch-case to handle different cabinet types
cabinetType = 2; % 1 for Closed Cabinet, 2 for Bass Reflex, 3 for Passive Radiator

switch cabinetType
    case 1
        select = 'Closed Cabinet';
        VB = 4.0e-3; % Cabinet volume (m^3)
        R = 1; % Monitoring distance (m)

        % Calculations specific to closed cabinet design
        PerformClosedCabinetCalculations;

    case 2
        select = 'Bass Reflex';
        VB = 4.0e-3; % Cabinet volume (m^3)
        R = 1; % Monitoring distance (m)

        % Calculations specific to bass reflex design
        PerformBassReflexCalculations;

    case 3
        select = 'Passive Radiator';
        VB = 4.0e-3; % Cabinet volume (m^3)
        R = 1; % Monitoring distance (m)

        % Calculations specific to passive radiator design
        PerformPassiveRadiatorCalculations;
end

% Plotting and annotations
title(select)
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB SPL)')
axis([10 1000 50 100])
grid on

% Closed Cabinet Calculations as a function
function [f, LT] = PerformClosedCabinetCalculations(f, s, UG, RE, SD, BL, rho, c, VB, MMS, CMS, RMS, pREF, R)
    FA = (BL*UG)/(RE*SD); % Acoustical force
    CAB = VB/(rho*c^2); % Box compliance
    RAE = BL^2/(RE*SD^2); % Electrical DC resistance equivalent
    MAS = MMS/(SD^2); % Driver moving mass
    CAS = CMS*SD^2; % Driver compliance
    RAS = RMS/SD^2; % Driver mechanical loss
    qF = FA./(RAE + s.*MAS + 1./(s.*CAS) + RAS + 1./(s.*CAB));
    pF = rho*s.*qF/(2*pi*R); % Volume velocity
    LT = 20*log10(abs(pF)/pREF); % Sound pressure level in dB SPL
end


% Bass Reflex Calculations as a function
function [f, LT, LD, LP] = PerformBassReflexCalculations(f, s, UG, RE, SD, BL, rho, c, VB, MMS, CMS, RMS, pREF, R)
    FA = (BL*UG)/(RE*SD); % Acoustical force
    CAB = VB/(rho*c^2); % Box compliance
    RAE = BL^2/(RE*SD^2); % Electrical DC resistance equivalent
    MAS = MMS/(SD^2); % Driver moving mass
    CAS = CMS*SD^2; % Driver compliance
    RAS = RMS/SD^2; % Driver mechanical loss
    RP = 0.015; % Port radius (m)
    SP = pi*RP^2; % Port area (m2)
    LX = 0.200; % Port length (m)
    MMP = (rho*SP)*(LX + 1.5*sqrt(SP/pi)); % Added mass of air due to the port
    MAP = MMP/SP^2; % Effective moving mass of the port
    qF = FA./(RAE + s.*MAS + 1./(s.*CAS) + RAS + 1./(s.*(CAB + 1./(s.*MAP))));
    qP = -qF.*(1./(s.*CAB))./(1./(s.*CAB) + s.*MAP);
    pT = rho*s.*(qF+qP)/(2*pi*R); % Total volume velocity
    pF = rho*s.*qF/(2*pi*R); % Driver volume velocity
    pP = rho*s.*qP/(2*pi*R); % Port volume velocity
    LT = 20*log10(abs(pT)/pREF); % Total sound pressure level
    LD = 20*log10(abs(pF)/pREF); % Driver sound pressure level
    LP = 20*log10(abs(pP)/pREF); % Port sound pressure level
end


% Passive Radiator Calculations as a function
function [f, LT, LF, LP] = PerformPassiveRadiatorCalculations(f, s, UG, RE, SD, BL, rho, c, VB, MMS, CMS, RMS, pREF, R)
    FA = (BL*UG)/(RE*SD); % Acoustical force
    CAB = VB/(rho*c^2); % Box compliance
    RAE = BL^2/(RE*SD^2); % Electrical DC resistance equivalent
    MAS = MMS/(SD^2); % Driver moving mass
    CAS = CMS*SD^2; % Driver compliance
    RAS = RMS/SD^2; % Driver mechanical loss
    MMP = 2*MMS; % Slave moving mass (for passive radiator)
    CMP = CMS; % Slave compliance (for passive radiator)
    RMP = RMS; % Slave mechanical loss (for passive radiator)
    MAP = MMP/SD^2; % Moving mass of the passive radiator
    CAP = CMP*SD^2; % Suspension compliance of the passive radiator
    RAP = RMP/SD^2; % Mechanical loss of the passive radiator
    qF = FA./(RAE + s.*MAS + 1./(s.*CAS) + RAS + 1./(s.*(CAB + 1./(s.*MAP + 1./(s.*CAP) + RAP))));
    qP = -qF.*(1

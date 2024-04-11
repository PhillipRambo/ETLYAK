udgangspunkt for simulerings milj

% Loudspeaker cabinet simulator
clear all
rho=1.18; % Air mass density (kg/m3).
c=345; % Speed of sound (m/s).
pREF=20e-6; % Reference sound pressure (Pa).
UG=2.83; % Voltage across terminals (V).
f=10:1000; % Frequency scale (Hz).
s=j*2*pi*f; % Laplace operator (rad/s).
% Loudspeaker parameters.
RE=5.6; % Coil DC resistance (ohm).
LE=0.65e-3; % Coil self-inductance (H).
BL=5.1; % Force factor (N/A).
MMS=0.007; % Diaphragm mass (kg).
CMS=1.45e-3; % Suspension compliance (m/N).
RMS=0.6; % Mechanical losses (Ns/m).
SD=0.0054; % Diaphragm area (m2).
switch 2
case {1}
select='Lukket kabinet';
VB=4.0e-3; % Cabinet volume (m3).
R=1; % Monitoring distance (m).
FA=(BL*UG/(RE*SD)); % Acoustical force.
CAB=VB/(rho*c^2); % Box compliance.
RAE=(BL)^2/(RE*SD^2); % Electrical DC resistance.
MAS=MMS/(SD^2); % Driver moving mass.
CAS=CMS*SD^2; % Driver compliance.
RAS=RMS/SD^2; % Driver mechanical loss.
RAL=10e4; % Cabinet losses (Ns/m5).
qF=FA./(RAE+s*MAS+1./(s*CAS)+RAS+1./(s*CAB));
pF=rho*s.*qF/(2*pi*R); % Volume velocity (m3/s).
LT=20*log10(abs(pF)/pREF); % Sound pressure (dB SPL).
semilogx(f,LT,'-r','LineWidth',2)
CMB=VB/(rho*c^2*SD^2);
CRES=CMS*CMB/(CMS+CMB);
FC=1/(2*pi*sqrt(MMS*CRES));
QTS=sqrt(MMS/CRES)/(BL^2/RE+RMS);
text(11,97, ['VB=' num2str(VB) ' liter'])
text(11,95, ['FC=' num2str(round(FC)) ' Hz'])
text(11,93, ['QTS=' num2str(round(100*QTS)/100)])
text(11,87, ['UG=' num2str(round(100*UG)/100) ' V'])
text(11,85, ['R=' num2str(round(100*R)/100) ' m'])
% -----------------------------------------------------
case {2}
select='Basreflex';
VB=4.0e-3; % Cabinet volume (m3).
R=1; % Monitoring distance (m).
FA=(BL*UG/(RE*SD)); % Acoustical force.
CAB=VB/(rho*c^2); % Box compliance.
RAE=(BL)^2/(RE*SD^2); % Electrical DC resistance.
MAS=MMS/(SD^2); % Driver moving mass (kg).
CAS=CMS*SD^2; % Driver compliance (m/N).
RAS=RMS/SD^2; % Driver mechanical loss (Ns/m).
RAL=10e4; % Cabinet losses (Ns/m5).
RP=0.015; % Port radius (m).
SP=pi*RP^2; % Port area (m2).
LX=0.200; % Port length (m).
MMP=(rho*SP)*(LX+1.5*sqrt(SP/pi));
MAP=MMP/SP^2; % Driver moving mass (kg).
qF=FA./(RAE+s*MAS+1./(s*CAS)+RAS+1./(s*CAB+1./(s*MAP)));
qP=-qF.*(1./(s*CAB))./(1./(s*CAB)+s*MAP);
pT=rho*s.*(qF+qP)/(2*pi*R); % Total volume velocity (m3/s).
pF=rho*s.*qF/(2*pi*R); % Driver volume velocity (m3/s).
pP=rho*s.*qP/(2*pi*R); % Port volume velocity (m3/s).
LT=20*log10(abs(pT)/pREF); % Total sound pressure (dB).
LD=20*log10(abs(pF)/pREF); % Driver sound pressure (dB).
LP=20*log10(abs(pP)/pREF); % Port sound pressure (dB).
semilogx(f,LT,'-r','LineWidth',2)
hold on
semilogx(f,LD,'-g','LineWidth',1)
semilogx(f,LP,'-b','LineWidth',1)
hold off
legend('Samlet','Enhed','Port')
CMB=VB/(rho*c^2*SD^2);
CRES=CMS*CMB/(CMS+CMB);
FC=1/(2*pi*sqrt(MMS*CRES));
FP=1/(2*pi*sqrt(MAP*CAB));
QTS=sqrt(MMS/CRES)/(BL^2/RE+RMS);
QP=RAE*sqrt(CAB/MAP);
text(11,97, ['VB=' num2str(1000*VB) ' liter'])
text(11,95, ['FC=' num2str(round(FC)) ' Hz'])
text(11,93, ['FP=' num2str(round(FP)) ' Hz'])
text(11,91, ['QTS=' num2str(round(100*QTS)/100)])
text(11,89, ['QP=' num2str(round(100*QP)/100)])
text(11,85, ['UG=' num2str(round(UG)) ' V'])
text(11,83, ['R=' num2str(round(R)) ' m'])
text(11,79, ['RP=' num2str(RP) ' m'])
text(11,77, ['LP=' num2str(LX) ' m'])
% -----------------------------------------------------
case {3}
select='Passiv slave';
VB=4.0e-3; % Cabinet volume (m3).
R=1; % Monitoring distance (m).
FA=(BL*UG/(RE*SD)); % Acoustical force.
CAB=VB/(rho*c^2); % Box compliance.
RAE=(BL)^2/(RE*SD^2); % Electrical DC resistance.
MMP=2*MMS; % Slave moving mass (kg).
CMP=CMS; % Slave compliance (m/N).
RMP=RMS; % Slave mechanical loss (Ns/m).
SP=SD; % Slave area (m2).
MAS=MMS/SP^2; % Driver moving mass.
CAS=CMS*SP^2; % Driver suspension compliance.
RAS=RMS/SP^2; % Driver mechanical loss.
MAP=MMP/SP^2; % Port moving mass.
CAP=CMP*SP^2; % Port suspension compliance.
RAP=RMP/SP^2; % Port mechanical loss.
qF=FA./(RAE+s*MAS+1./(s*CAS)+RAS+1./(s*CAB+1./(s*MAP+1./(s*CAP)+RAP)));

qP=-qF.*(1./(s*CAB))./(1./(s*CAB)+s*MAP+1./(s*CAP)+RAP);
pT=rho*s.*(qF+qP)/(2*pi*R); % Total volume velocity (m3/s).
pF=rho*s.*qF/(2*pi*R); % Driver volume velocity (m3/s).
pP=rho*s.*qP/(2*pi*R); % Port volume velocity (m3/s).
LT=20*log10(abs(pT)/pREF); % Total sound pressure (dB).
LF=20*log10(abs(pF)/pREF); % Driver sound pressure (dB).
LP=20*log10(abs(pP)/pREF); % Slave sound pressure (dB).
% semilogx(f,LT,'-r','Linewidth',2, f,LF,'-g', f,LP,'-b')
semilogx(f,LT,'-r','Linewidth',2)
hold on
semilogx(f,LF,'-g')
semilogx(f,LP,'-b')
hold off
legend('Samlet','Enhed','Slave')
text(11,97, ['VB=' num2str(1000*VB) ' liter'])
CMB=VB/(rho*c^2*SD^2);
CRES=CMS*CMB/(CMS+CMB);
FC=1/(2*pi*sqrt(MMS*CRES));
FP=1/(2*pi*sqrt(MAP*CAB));
QTS=sqrt(MMS/CRES)/(BL^2/RE+RMS);
QP=(1/RAP)*sqrt(MAP/CAP);
text(11,95, ['FC=' num2str(round(FC)) ' Hz'])
text(11,93, ['FP=' num2str(round(FP)) ' Hz'])
text(11,91, ['QTS=' num2str(round(100*QTS)/100)])
text(11,89, ['QP=' num2str(round(100*QP)/100)])
text(11,85, ['UG=' num2str(round(100*UG)/100) ' V'])
text(11,83, ['R=' num2str(round(100*R)/100) ' m'])
end
title(select)
xlabel('Frekvens (Hz)')
ylabel('Amplitude (dB SPL)')
axis([10 1000 50 100] )
grid on


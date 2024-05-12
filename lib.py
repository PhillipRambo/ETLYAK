from dataclasses import dataclass
import numpy as np
import numpy as np


@dataclass
class TsParams:
    bl: float # Force Factor in Tesla-meters
    mms: float  # Moving Mass in grams
    cms: float  # Suspension Compliance in mm/N
    Qes: float  # Electrical Q, dimensionless
    Qms: float  # Mechanical Q, dimensionless
    Re: float  # Voice Coil DC Resistance in ohms
    rms: float  # Mechanical Losses in Ns/m
    cms: float # Suspension compliance af slaven.
    rms: float # Mekaniske tab af slaven.
    sd: float # Slavens areal.




@dataclass
class SpeakerUnit:
    params: TsParams
    surface_area: float


def speaker_from_blue_planet_parameters(impedance: float, xmax: float, fres: float, bl: float, Le: float, Re: float, Qms: float, Qes: float, Qts: float, Vas:float, Sd: float, Mms: float, Cms: float, Rms: float):
    ts_params_local = TsParams(bl, Mms, Cms, Qes, Qms, Re, Rms, Sd)
    surface_area = Sd
    return SpeakerUnit(ts_params_local, surface_area)
        

@dataclass
class Cabinet:
    width: float
    height: float
    depth: float
    volume: float
    port_size: float
    port_N: int
    port_length : float

@dataclass
class BassReflex:
    unit: SpeakerUnit
    cabinet: Cabinet

@dataclass
class PassiveSlave:
    unit: SpeakerUnit
    cabinet: Cabinet


    
"""
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
"""

def simulate_bass_reflex(bassreflex: BassReflex, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 345  # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    sd = bassreflex.unit.surface_area

    ts = bassreflex.unit.params

    fa = (ts.Qes * ts.Qms) / (ts.Re * sd) # Acoustical force
    cab = bassreflex.cabinet.volume / (rho * c**2) # Box compliance
    rae = ts.bl**2 / (ts.Re * sd**2) # Electrical DC resistance equivalent
    mas = ts.mms / (sd**2) # Driver moving mass
    cas = ts.cms * sd**2 # Driver compliance
    ras = ts.rms / sd**2 # Driver mechanical loss
    rp = bassreflex.cabinet.port_size # Port radius (m)
    sp = np.pi * rp**2 # Port area (m2)
    lx = bassreflex.cabinet.port_length # Port length (m)
    mmp = (rho * sp) * (lx + 1.5 * np.sqrt(sp / np.pi)) # Added mass of air due to the port
    map = mmp / sp**2 # Effective moving mass of the port
    qf = fa / (rae + s * mas + 1 / (s * cas) + ras + 1 / (s * (cab + 1 / (s * map))))
    qp = -qf * (1 / (s * cab)) / (1 / (s * cab) + s * map) # Port volume velocity
    pt = rho * s * (qf + qp) / (2 * np.pi * ts.Re) # Total volume velocity
    pf = rho * s * qf / (2 * np.pi * ts.Re) # Driver volume velocity
    pp = rho * s * qp / (2 * np.pi * ts.Re) # Port volume velocity
    lt = 20 * np.log10(np.abs(pt) / pREF) # Total sound pressure level
    ld = 20 * np.log10(np.abs(pf) / pREF) # Driver sound pressure level
    lp = 20 * np.log10(np.abs(pp) / pREF) # Port sound pressure level
    return f, lt, ld, lp



def passiveslave_simulation(passive_slave: PassiveSlave, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 345     # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    ts = passive_slave.unit.params
    sd = passive_slave.unit.surface_area

    fa = (ts.Qes * ts.Qms) / (ts.Re * sd) 
    cab = passive_slave.cabinet.volume / (rho * c**2)  # Box compliance for the slave
    rae = ts.bl**2 / (ts.Re * sd**2)  # Electrical DC resistance equivalent

    mmp = ts.mmp
    cmp = ts.cms
    rmp = ts.rms
    sp = ts.sd

    mas = ts.mms / (sp**2)
    cas = ts.cms * sp**2
    ras = ts.rms / sp**2
    map = mmp / sp**2
    cap = cmp * sp**2
    rap = rmp / sp**2

    qF = fa / (rae + s * mas + 1 / (s * cas) + ras + 1 / (s * cab + 1 / (s * map + 1 / (s * cap) + rap)))
    qP = -qF * (1 / (s * cab)) / (1 / (s * cab) + s * map + 1 / (s * cap) + rap)

    pT = rho * s * (qF + qP) / (2 * np.pi * ts.Re)
    pF = rho * s * qF / (2 * np.pi * ts.Re)
    pP = rho * s * qP / (2 * np.pi * ts.Re)

    lt = 20 * np.log10(np.abs(pT) / pREF)
    lf = 20 * np.log10(np.abs(pF) / pREF)
    lp = 20 * np.log10(np.abs(pP) / pREF)

    return f, lt, lf, lp





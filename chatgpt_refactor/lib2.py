from dataclasses import dataclass
import numpy as np
import numpy as np
import math
from abc import ABC
from enum import Enum

class SpeakerType(ABC):
    pass

class SimulationType(Enum):
    BASS_REFLEX = 1
    PASSIVE_SLAVE = 2
    BANDPASS_6TH_ORDER = 3

@dataclass
class TsParams:
    Bl: float # Force Factor in Tesla-meters
    Mms: float  # Moving Mass in grams
    Cms: float  # Suspension Compliance in mm/N
    Qes: float  # Electrical Q, dimensionless
    Qms: float  # Mechanical Q, dimensionless
    Re: float  # Voice Coil DC Resistance in ohms
    Rms: float  # Mechanical Losses in Ns/m
    Cms: float # Suspension compliance
    Rms: float # Mekaniske tab
    Sd: float # Membrane area


@dataclass
class SpeakerUnit:
    params: TsParams


def unit_from_blue_planet_parameters(impedance: float, xmax: float, fres: float, bl: float, Le: float, Re: float, Qms: float, Qes: float, Qts: float, Vas:float, Sd: float, Mms: float, Cms: float, Rms: float):
    ts_params_local = TsParams(bl, Mms, Cms, Qes, Qms, Re, Rms, Sd)
    return SpeakerUnit(ts_params_local)
        

@dataclass
class Cabinet:
    width: float
    height: float
    depth: float
    volume: float

    def __init__(self, width, height, depth):
        self.width = width
        self.height = height
        self.depth = depth
        self.volume = width * height * depth
    
    def __init__(self, volume):
        self.volume = volume # Volume in liters
        self.width = math.pow(volume, 1/3) # Assume cubic shape
        self.height = math.pow(volume, 1/3)
        self.depth = math.pow(volume, 1/3)

@dataclass
class Port:
    radius: float # Port radius in cm
    length: float # Port length in cm

    def __init__(self, radius, length):
        self.radius = radius
        self.length = length
        self.area = math.pi * radius**2

@dataclass
class PassiveUnit:
    Cas: float
    Mas: float
    Ras: float

@dataclass
class BassReflex(SpeakerType):
    unit: SpeakerUnit
    cabinet: Cabinet
    port: Port

    def __init__(self, unit: SpeakerUnit, cabinet: Cabinet = Cabinet(8), port: Port = Port(2, 10)):
        self.unit = unit
        self.cabinet = Cabinet(8)
        self.port = Port(2, 10)

    def __str__(self) -> str:
        return "Bass Reflex"

@dataclass
class PassiveSlave(SpeakerType):
    unit: SpeakerUnit
    slave: PassiveUnit
    cabinet: Cabinet

    def __init__(self, unit: SpeakerUnit):
        self.unit = unit
        self.cabinet = Cabinet(8)
        self.slave = PassiveUnit(self.unit.params.Cms, self.unit.params.Mms, self.unit.params.Rms)

    def __str__(self) -> str:
        return "Passive Slave"


@dataclass
class Bandpass6thOrder(SpeakerType):
    unit: SpeakerUnit
    front_cabinet: Cabinet
    back_cabinet: Cabinet
    front_port: Port
    back_ports: Port

    def __init__(self, unit: SpeakerUnit):
        self.unit = unit
        self.front_cabinet = Cabinet(1)
        self.back_cabinet = Cabinet(30)
        self.front_port = Port(3, 10)
        self.back_ports = Port(3, 17)
    
    def __str__(self) -> str:
        return "6th Order Bandpass"


def simulate_bass_reflex(bassreflex: BassReflex, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 343  # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    ts = bassreflex.unit.params
    print(ts)

    Ug = 1 # Amplitude of the input signal
    Ug_eq = (ts.Bl)/(ts.Re*ts.Sd) * Ug # Equivalent input signal for acoustical circuit
    r = 1 # Distance from the speaker to the listener

    #* Driver acoustical impedance
    Zrae = (ts.Bl**2) / (ts.Re * ts.Sd**2) # Electrical DC resistance equivalent
    Zmas = s * (ts.Mms) / (ts.Sd**2) # Mechanical impedance of the driver
    Zcas = 1 / (s * ts.Cms * ts.Sd**2) # Compliance impedance of the driver
    Zras = (ts.Rms) / (ts.Sd**2) # Mechanical impedance of the driver

    #* Front acoustical load impedance
    Zaf = 0 # unit presents no load
    #* Rear acoustical load impedance
    Vbr = bassreflex.cabinet.volume * 1e-3 # Volume of the rear chamber in m^3
    Car = Vbr / (rho * c**2) # Compliance of the rear chamber
    Zcar = 1 / (s * Car) # Acoustical impedance of the rear chamber
    Spr = np.pi * (bassreflex.port.radius*1e-2)**2 # Area of the rear ports
    Lpr = bassreflex.port.length*1e-2 # Length of the rear ports
    Mapr = rho/Spr * (Lpr + 1.5 * np.sqrt(Spr / np.pi)) # Added mass of air due to the rear ports
    Zmar = s * Mapr # Mechanical impedance of the rear ports

    Zab = (Zcar * Zmar) / (Zcar + Zmar) # Total acoustical impedance of the back chamber

    #* Air flow through circuit
    q = (Ug_eq)/(Zrae + Zmas + Zcas + Zras + Zaf + Zab)

    #* Sound pressure level
    p_factor = s*rho/(2 * np.pi * r)
    # current divider between Zcaf and Zmaf - only current in Zmaf is converted to sound pressure
    pf = p_factor * q
    # current divider between Zcar and Zmar - only current in Zmar is converted to sound pressure
    pr = p_factor * (-q) * Zcar / (Zcar + Zmar)

    p_total = pf + pr

    splF = 20 * np.log10(np.abs(pf) / pREF)
    splR = 20 * np.log10(np.abs(pr) / pREF)
    splT = 20 * np.log10(np.abs(p_total) / pREF)

    return f, splT, splF, splR


def simulate_passive_slave(passive_slave: PassiveSlave, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 343  # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    ts = passive_slave.unit.params
    print(ts)

    Ug = 1 # Amplitude of the input signal
    Ug_eq = (ts.Bl)/(ts.Re*ts.Sd) * Ug # Equivalent input signal for acoustical circuit
    r = 1 # Distance from the speaker to the listener

    #* Driver acoustical impedance
    Zrae = (ts.Bl**2) / (ts.Re * ts.Sd**2) # Electrical DC resistance equivalent
    Zmas = s * (ts.Mms) / (ts.Sd**2) # Mechanical impedance of the driver
    Zcas = 1 / (s * ts.Cms * ts.Sd**2) # Compliance impedance of the driver
    Zras = (ts.Rms) / (ts.Sd**2) # Mechanical impedance of the driver

    #* Front acoustical load impedance
    Zaf = 0 # unit presents no load
    #* Rear acoustical load impedance
    Vbr = passive_slave.cabinet.volume * 1e-3 # Volume of the rear chamber in m^3
    Car = Vbr / (rho * c**2) # Compliance of the rear chamber
    Zcar = 1 / (s * Car) # Acoustical impedance of the rear chamber
    Zslave = passive_slave.slave.Ras + s * passive_slave.slave.Mas + 1 / (s * passive_slave.slave.Cas) # Acoustical impedance of the slave

    Zab = (Zcar * Zslave) / (Zcar + Zslave) # Total acoustical impedance of the back chamber

    #* Air flow through circuit
    q = (Ug_eq)/(Zrae + Zmas + Zcas + Zras + Zaf + Zab)

    #* Sound pressure level
    p_factor = s*rho/(2 * np.pi * r)
    # current divider between Zcaf and Zmaf - only current in Zmaf is converted to sound pressure
    pf = p_factor * q
    # current divider between Zcar and Zmar - only current in Zmar is converted to sound pressure
    pr = p_factor * (-q) * Zcar / (Zcar + Zslave)

    p_total = pf + pr

    splF = 20 * np.log10(np.abs(pf) / pREF)
    splR = 20 * np.log10(np.abs(pr) / pREF)
    splT = 20 * np.log10(np.abs(p_total) / pREF)

    return f, splT, splF, splR


def simulate_6thorderbandpass(bandpass: Bandpass6thOrder, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 343  # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    ts = bandpass.unit.params
    print(ts)

    Ug = 1 # Amplitude of the input signal
    Ug_eq = (ts.Bl)/(ts.Re*ts.Sd) * Ug # Equivalent input signal for acoustical circuit
    r = 1 # Distance from the speaker to the listener

    #* Driver acoustical impedance
    Zrae = (ts.Bl**2) / (ts.Re * ts.Sd**2) # Electrical DC resistance equivalent
    Zmas = s * (ts.Mms) / (ts.Sd**2) # Mechanical impedance of the driver
    Zcas = 1 / (s * ts.Cms * ts.Sd**2) # Compliance impedance of the driver
    Zras = (ts.Rms) / (ts.Sd**2) # Mechanical impedance of the driver

    #* Front acoustical load impedance
    Vbf = bandpass.front_cabinet.volume * 1e-3 # Volume of the front chamber in m^3
    Caf = Vbf / (rho * c**2) # Compliance of the front chamber
    Zcaf = 1 / (s * Caf) # Acoustical impedance of the front chamber
    Spf = np.pi * (bandpass.front_port.radius*1e-2)**2 # Area of the front port
    Lpf = bandpass.front_port.length*1e-2 # Length of the front port
    Mapf = rho/Spf * (Lpf + 1.5 * np.sqrt(Spf / np.pi)) # Added mass of air due to the front port
    Zmaf = s * Mapf # Mechanical impedance of the front port

    Zaf = (Zcaf * Zmaf) / (Zcaf + Zmaf) # Total acoustical impedance of the front chamber

    #* Rear acoustical load impedance
    Vbr = bandpass.back_cabinet.volume * 1e-3 # Volume of the rear chamber in m^3
    Car = Vbr / (rho * c**2) # Compliance of the rear chamber
    Zcar = 1 / (s * Car) # Acoustical impedance of the rear chamber
    Spr = np.pi * (bandpass.back_ports.radius*1e-2)**2 # Area of the rear ports
    Lpr = bandpass.back_ports.length*1e-2 # Length of the rear ports
    Mapr = rho/Spr * (Lpr + 1.5 * np.sqrt(Spr / np.pi)) # Added mass of air due to the rear ports
    Zmar = s * Mapr # Mechanical impedance of the rear ports

    Zab = (Zcar * Zmar) / (Zcar + Zmar) # Total acoustical impedance of the back chamber

    #* Air flow through circuit
    q = (Ug_eq)/(Zrae + Zmas + Zcas + Zras + Zaf + Zab)

    #* Sound pressure level
    p_factor = s*rho/(2 * np.pi * r)
    # current divider between Zcaf and Zmaf - only current in Zmaf is converted to sound pressure
    pf = p_factor * q * Zcaf / (Zcaf + Zmaf)
    # current divider between Zcar and Zmar - only current in Zmar is converted to sound pressure
    pr = p_factor * (-q) * Zcar / (Zcar + Zmar)

    p_total = pf + pr

    splF = 20 * np.log10(np.abs(pf) / pREF)
    splR = 20 * np.log10(np.abs(pr) / pREF)
    splT = 20 * np.log10(np.abs(p_total) / pREF)

    return f, splT, splF, splR

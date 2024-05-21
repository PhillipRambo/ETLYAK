import numpy as np
import matplotlib.pyplot as plt

def passiveradiator_6thorderbandpass_simulation(Slave, bandpass, frequency_range=(10, 1000)):
    rho = 1.18  # Air mass density (kg/m^3)
    c = 343  # Speed of sound (m/s)
    pREF = 20e-6  # Reference sound pressure (Pa)

    f = np.arange(frequency_range[0], frequency_range[1] + 1)
    s = 1j * 2 * np.pi * f

    # Access parameters directly from the unit instances
    ts1 = bandpass.unit1
    ts2 = bandpass.unit2

    Ug = 1  # Amplitude of the input signal
    Ug_eq = (ts1.Sd * Ug + ts2.Sd * Ug) / (ts1.Bl * ts1.Re + ts2.Bl * ts2.Re)  # Equivalent input signal for acoustical circuit
    r = 1  # Distance from the speaker to the listener

    # Driver acoustical impedance
    Zrae1 = (ts1.Bl**2) / (ts1.Re * ts1.Sd**2)
    Zmas1 = s * (ts1.Mms) / (ts1.Sd**2)
    Zcas1 = 1 / (s * ts1.Cms * ts1.Sd**2)
    Zras1 = (ts1.Rms) / (ts1.Sd**2)

    Zrae2 = (ts2.Bl**2) / (ts2.Re * ts2.Sd**2)
    Zmas2 = s * (ts2.Mms) / (ts2.Sd**2)
    Zcas2 = 1 / (s * ts2.Cms * ts2.Sd**2)
    Zras2 = (ts2.Rms) / (ts2.Sd**2)

    # Front acoustical load impedance
    Vbf = bandpass.front_cabinet.volume * 1e-3  # Volume of the front chamber in m^3
    Cab = Vbf / (rho * c**2)  # Compliance of the front chamber
    Map = Slave.Mmp / Slave.Sp**2
    Cap = Slave.Cmp / Slave.Sp**2
    Rap = Slave.Rmp / Slave.Sp**2
    Zcaf = 1 / ((s * Cab) + 1 / (s * Map + 1 / (s * Cap) + Rap))

    # Rear acoustical load impedance
    Vbr = bandpass.back_cabinet.volume * 1e-3  # Volume of the rear chamber in m^3
    Car = Vbr / (rho * c**2)  # Compliance of the rear chamber
    Zcar = 1 / (s * Car)  # Acoustical impedance of the rear chamber

    # Bass reflex ports
    Spr1 = np.pi * (bandpass.back_ports[0].radius * 1e-2)**2  # Area of the first rear port
    Lpr1 = bandpass.back_ports[0].length * 1e-2  # Length of the first rear port
    Mapr1 = rho / Spr1 * (Lpr1 + 1.5 * np.sqrt(Spr1 / np.pi))  # Added mass of air due to the first rear port
    Zmar1 = 1 / 2 * s * Mapr1  # Mechanical impedance of the first rear port

    Spr2 = np.pi * (bandpass.back_ports[1].radius * 1e-2)**2  # Area of the second rear port
    Lpr2 = bandpass.back_ports[1].length * 1e-2  # Length of the second rear port
    Mapr2 = rho / Spr2 * (Lpr2 + 1.5 * np.sqrt(Spr2 / np.pi))  # Added mass of air due to the second rear port
    Zmar2 = 1 / 2 * s * Mapr2  # Mechanical impedance of the second rear port

    # Total acoustical impedance of the back chamber
    Zabr1 = (Zcar * Zmar1) / (Zcar + Zmar1)
    Zabr2 = (Zcar * Zmar2) / (Zcar + Zmar2)
    Zab = Zabr1 + Zabr2

    # Airflow through the whole circuit
    qF = Ug_eq / (Zrae1 + Zmas1 + Zcas1 + Zras1 + Zrae2 + Zmas2 + Zcas2 + Zras2 + Zcaf + Zab)

    # Airflow through the passive radiator
    qP = - (1 / (s * Cab)) / (1 / (s * Cab) + s * Map + 1 / (s * Cap) + Rap) * qF

    # Volume velocities
    pT = rho * s * (qF + qP) / (2 * np.pi * r)  # Total volume velocity
    pF = rho * s * qF / (2 * np.pi * r)  # Driver volume velocity
    pP = rho * s * qP / (2 * np.pi * r)  # Passive radiator volume velocity

    # Sound pressure levels
    splF = 20 * np.log10(np.abs(pF) / pREF)  # Driver sound pressure level
    splP = 20 * np.log10(np.abs(pP) / pREF)  # Passive radiator sound pressure level
    splT = 20 * np.log10(np.abs(pT) / pREF)  # Total sound pressure level

    return f, splT, splF, splP

# Example usage (replace with actual parameters)
class PassiveUnit:
    def __init__(self, Mmp, Sp, Cmp, Rmp):
        self.Mmp = Mmp
        self.Sp = Sp
        self.Cmp = Cmp
        self.Rmp = Rmp

class Port:
    def __init__(self, radius, length):
        self.radius = radius
        self.length = length

class Cabinet:
    def __init__(self, volume):
        self.volume = volume

class UnitParams:
    def __init__(self, Sd, Bl, Re, Mms, Cms, Rms):
        self.Sd = Sd
        self.Bl = Bl
        self.Re = Re
        self.Mms = Mms
        self.Cms = Cms
        self.Rms = Rms

class Bandpass6thOrderPassiveSlave:
    def __init__(self, unit1, unit2, front_cabinet, back_cabinet, back_ports):
        self.unit1 = unit1
        self.unit2 = unit2
        self.front_cabinet = front_cabinet
        self.back_cabinet = back_cabinet
        self.back_ports = back_ports

# Define parameters
slave = PassiveUnit(Mmp=0.2, Sp=0.05, Cmp=0.001, Rmp=0.1)
unit1 = UnitParams(Sd=0.05, Bl=10, Re=6, Mms=0.02, Cms=1e-4, Rms=0.5)
unit2 = UnitParams(Sd=0.05, Bl=10, Re=6, Mms=0.02, Cms=1e-4, Rms=0.5)
front_cabinet = Cabinet(volume=20)
back_cabinet = Cabinet(volume=40)
back_ports = [Port(radius=0.05, length=0.1), Port(radius=0.05, length=0.1)]

bandpass = Bandpass6thOrderPassiveSlave(unit1, unit2, front_cabinet, back_cabinet, back_ports)

# Run simulation
frequency, spl_total, spl_driver, spl_radiator = passiveradiator_6thorderbandpass_simulation(slave, bandpass)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(frequency, spl_total, label='Total SPL')
plt.plot(frequency, spl_driver, label='Driver SPL')
plt.plot(frequency, spl_radiator, label='Radiator SPL')
plt.xscale('log')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Sound Pressure Level (dB)')
plt.title('Sound Pressure Levels')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()

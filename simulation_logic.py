import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lib import *

cabinet = Cabinet(width=20, height=30, depth=20, volume=0.02, port_size=10, port_N=1, port_length=10)
speaker_unit = speaker_from_blue_planet_parameters(impedance=4, xmax=8, fres=45, bl=8.7, Le=1.18, Re=3.5, Qms=4.37, Qes=0.49, Qts=0.44, Vas=6.8, Sd=127, Mms=38.64, Cms=0.3, Rms=2.57)
basreflex = BassReflex(speaker_unit, cabinet)


def plot_speaker_response():
    # No need for VB as a parameter if you're accessing it directly from the basreflex object
    frequency, total_sound_pressure, driver_pressure_level, port_pressure_level = simulate_bass_reflex(basreflex)
    ax.clear()
    ax.plot(frequency, total_sound_pressure, '-r', linewidth=2, label='Total SPL')
    ax.plot(frequency, driver_pressure_level, '-g', linewidth=2, label='Driver SPL')
    ax.plot(frequency, port_pressure_level, '-b', linewidth=2, label='Port SPL')
    ax.legend(loc='best')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude (dB SPL)')
    ax.set_title(f'Loudspeaker Simulation (Volume {basreflex.cabinet.volume} L)')
    ax.grid(True)
    canvas.draw()

def on_volume_entry_change(*args):
    # Get the volume from the input box
    volume_str = volume_entry.get()
    if volume_str:
        try:
            # Try to convert it to a float
            VB = float(volume_str)
            # Update the slider's position
            volume_slider.set(VB)
            # Update the cabinet's volume
            basreflex.cabinet.volume = VB
            # Redraw the plot
            plot_speaker_response()
        except ValueError:
            # Handle the exception if the input is not a valid float
            print("Please enter a valid number for volume.")


def on_slider_move(event):
    VB = volume_slider.get()
    basreflex.cabinet.volume = VB
    volume_entry.delete
    volume_entry.insert(0, str(VB))
    plot_speaker_response()
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lib import *

cabinet = Cabinet(width=20, height=30, depth=20, volume=0.02, port_size=10, port_N=1, port_length=10)
speaker_unit = speaker_from_blue_planet_parameters(impedance=4, xmax=8, fres=45, bl=8.7, Le=1.18, Re=3.5, Qms=4.37, Qes=0.49, Qts=0.44, Vas=6.8, Sd=127, Mms=38.64, Cms=0.3, Rms=2.57)
basreflex = BassReflex(speaker_unit, cabinet)


def plot_speaker_response(VB):
    frequency, total_sound_pressure, driver_pressure_level, port_pressure_level = simulate_bass_reflex(basreflex)
    ax.clear()
    ax.plot(frequency, total_sound_pressure, '-r', linewidth=2)
    ax.plot(frequency, driver_pressure_level, '-g', linewidth=2)
    ax.plot(frequency, port_pressure_level, '-b', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude (dB SPL)')
    ax.set_title(f'Loudspeaker Simulation (Volume {VB} L)')
    #print default cabinet values
    ax.text(0.5, 0.5, f'Width: {cabinet.width} cm\nHeight: {cabinet.height} cm\nDepth: {cabinet.depth} cm\nVolume: {cabinet.volume} L\nPort Size: {cabinet.port_size} cm^2\nPort N: {cabinet.port_N}\nPort Length: {cabinet.port_length} cm', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=10)
    ax.grid(True)
    canvas.draw()

def on_slider_move(event):
    VB = volume_entry.get()
    if VB.strip():  # Check if the entry is not empty
        VB = int(VB)
        plot_speaker_response(VB)

def submit():
    # Get the values from entry boxes
    width_str = width_entry.get()
    height_str = height_entry.get()
    depth_str = depth_entry.get()
    volume_str = volume_entry.get()
    port_size_str = port_size_entry.get()
    port_N_str = port_N_entry.get()
    port_length_str = port_length_entry.get()

    # Validate the input values
    if width_str.strip() and height_str.strip() and depth_str.strip() and volume_str.strip() and port_size_str.strip() and port_N_str.strip() and port_length_str.strip():
        # Convert input values to float
        width = float(width_str)
        height = float(height_str)
        depth = float(depth_str)
        volume = float(volume_str)
        port_size = float(port_size_str)
        port_N = int(port_N_str)
        port_length = float(port_length_str)

        # Create the cabinet object
        cabinet = Cabinet(width=width, height=height, depth=depth, volume=volume, port_size=port_size, port_N=port_N, port_length=port_length)

        print(cabinet)
        update_speaker()
    else:
        # Display error message in red
        #display the error message on the gui
        error_message = tk.Label(input_frame, text="Error: Please enter all values", font=('Helvetica', 12, 'bold'), background='#f0f0f0', fg='red')
        error_message.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky='ew')
        error_message.after(3000, error_message.destroy)
        
        print("Error: Please enter all values")

def update_speaker():
    # Clear existing speaker
    speaker_canvas.delete("all")

    # Get dimensions from entry boxes
    width = float(width_entry.get())
    height = float(height_entry.get())
    depth = float(depth_entry.get())

    # Draw simple speaker illustration
    x0 = 10
    y0 = 10
    x1 = x0 + width
    y1 = y0 + height
    x2 = x0 + width / 2
    y2 = y1 + depth
    speaker_canvas.create_rectangle(x0, y0, x1, y1, outline="black", fill="grey")
    speaker_canvas.create_polygon(x0, y1, x1, y1, x2, y2, outline="black", fill="grey")

# Create main window
window = tk.Tk()
window.title("Loudspeaker Simulation")
window.geometry("800x600")
window.configure(bg='#f0f0f0')

# Create matplotlib figure and canvas
fig = Figure(figsize=(6, 6), dpi=100)
ax = fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=window)
canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

# Create input frame
input_frame = tk.Frame(window, bg='#f0f0f0')
input_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=20, pady=20)

# Add input labels and entry boxes
volume_label = ttk.Label(input_frame, text="Volume (L):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
volume_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')
volume_entry = ttk.Entry(input_frame)
volume_entry.grid(row=0, column=1, padx=5, pady=5, sticky='ew')

width_label = ttk.Label(input_frame, text="Width (cm):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
width_label.grid(row=1, column=0, padx=5, pady=5, sticky='w')
width_entry = ttk.Entry(input_frame)
width_entry.grid(row=1, column=1, padx=5, pady=5, sticky='ew')

height_label = ttk.Label(input_frame, text="Height (cm):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
height_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
height_entry = ttk.Entry(input_frame)
height_entry.grid(row=2, column=1, padx=5, pady=5, sticky='ew')

depth_label = ttk.Label(input_frame, text="Depth (cm):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
depth_label.grid(row=3, column=0, padx=5, pady=5, sticky='w')
depth_entry = ttk.Entry(input_frame)
depth_entry.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

port_size_label = ttk.Label(input_frame, text="Port Size (cm^2):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
port_size_label.grid(row=4, column=0, padx=5, pady=5, sticky='w')
port_size_entry = ttk.Entry(input_frame)
port_size_entry.grid(row=4, column=1, padx=5, pady=5, sticky='ew')

port_N_label = ttk.Label(input_frame, text="Port N:", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
port_N_label.grid(row=5, column=0, padx=5, pady=5, sticky='w')
port_N_entry = ttk.Entry(input_frame)
port_N_entry.grid(row=5, column=1, padx=5, pady=5, sticky='ew')

port_length_label = ttk.Label(input_frame, text="Port Length (cm):", font=('Helvetica', 12, 'bold'), background='#f0f0f0')
port_length_label.grid(row=6, column=0, padx=5, pady=5, sticky='w')
port_length_entry = ttk.Entry(input_frame)
port_length_entry.grid(row=6, column=1, padx=5, pady=5, sticky='ew')

submit_button = ttk.Button(input_frame, text="Submit", command=submit)
submit_button.grid(row=7, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

# Create slider frame
slider_frame = tk.Frame(window, bg='#f0f0f0')
slider_frame.pack(side=tk.RIGHT, fill=tk.BOTH, padx=20, pady=20)

# Add slider for volume
volume_slider_label = tk.Label(slider_frame, text="Volume (L):", bg='#f0f0f0', fg='black', font=('Helvetica', 12, 'bold'))
volume_slider_label.pack(side=tk.TOP, padx=10, pady=10)

volume_slider = tk.Scale(slider_frame, from_=0, to=100, orient=tk.HORIZONTAL, length=300, bg='#f0f0f0', fg='blue', highlightbackground='#f0f0f0', troughcolor='#c0c0c0', command=on_slider_move)
volume_slider.pack(side=tk.BOTTOM, padx=10, pady=10)

# Create canvas for speaker illustration
speaker_canvas = tk.Canvas(slider_frame, width=150, height=200, bg='#f0f0f0', bd=0, highlightthickness=0)
speaker_canvas.pack()

window.mainloop()










''''
mport tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from lib import *

# Your simulate_loudspeaker function should be defined here

window = tk.Tk()
window.title("Graph Plotter")
window.geometry("800x600")
window.config(bg="white")

def plot(event): 
    VB = plot_slider.get() * 1e-3  # Convert from liters to cubic meters
    frequencies, sound_pressure = simulate_loudspeaker(VB=VB)

    fig = Figure(figsize=(5, 5), dpi=100)
    plot1 = fig.add_subplot(111)
    plot1.plot(frequencies, sound_pressure, '-r', linewidth=2)
    plot1.set_xlabel('Frequency (Hz)')
    plot1.set_ylabel('Amplitude (dB SPL)')
    plot1.set_title('Loudspeaker Simulation')
    plot1.grid(True)

    # Remove any existing canvas before drawing a new one
    for widget in window.winfo_children():
        if isinstance(widget, FigureCanvasTkAgg):
            widget.get_tk_widget().destroy()

    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    toolbar = NavigationToolbar2Tk(canvas, window)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    

    




plot_slider = tk.Scale(master=window, from_=1, to=100, orient=tk.HORIZONTAL, label="Cabinet Volume (L)", sliderlength=15, length=300) # slider for cabinet volume

plot_slider.bind("<ButtonRelease-1>", plot) # bind the plot function to the slider release event


## reconstruct plot after slider is released
plot_slider.pack(side=tk.BOTTOM, fill=tk.X)





tk.mainloop()

'''''


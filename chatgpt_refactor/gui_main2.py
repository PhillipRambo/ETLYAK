import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lib2 import unit_from_blue_planet_parameters, SpeakerType, BassReflex, Cabinet, Port, PassiveUnit, PassiveSlave, SimulationType, Bandpass6thOrder

UNIT = unit_from_blue_planet_parameters(impedance=4, xmax=8, fres=45, bl=8.7, Le=1.18, Re=3.5, Qms=4.37, Qes=0.49, Qts=0.44, Vas=6.8, Sd=127e-3, Mms=38.64e-3, Cms=0.3e-3, Rms=2.57)

def setup_gui():
    # Main window setup for input controls
    window = tk.Tk()
    window.title("Loudspeaker Parameter Input")
    window.geometry("300x400")
    window.configure(bg='#f0f0f0')

    # Separate window for plotting
    plot_window = tk.Toplevel(window)
    plot_window.title("Loudspeaker Simulation Plot")
    plot_window.geometry("800x600")

    speaker = BassReflex(UNIT, Cabinet(8), Port(2, 10))

    fig, ax, canvas = setup_plot_window(plot_window)

    setup_input_controls(window, speaker, ax, canvas)

    window.mainloop()

def setup_plot_window(window):
    fig = Figure(figsize=(8, 6), dpi=100)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    return fig, ax, canvas

def setup_input_controls(window, speaker: SpeakerType, ax, canvas):
    frame = ttk.Frame(window, padding="3 3 12 12")
    frame.pack(fill=tk.BOTH, expand=True)

    simulation_types = ['Bass Reflex', 'Passive Slave', '6th Order Bandpass']
    selected_type = tk.StringVar(value=simulation_types[0])

    ttk.Label(frame, text="Select Simulation Type:").pack()
    dropdown = ttk.Combobox(frame, textvariable=selected_type, values=simulation_types, state="readonly")
    dropdown.pack()

    params_frame = ttk.Frame(frame)
    params_frame.pack(fill=tk.BOTH, expand=True, pady=20)

    submit_button = ttk.Button(frame, text="Submit", command=lambda: submit(params_frame, speaker, ax, canvas))
    submit_button.pack()

    dropdown.bind("<<ComboboxSelected>>", lambda event: update_inputs(params_frame, selected_type.get(), speaker))
    update_inputs(params_frame, selected_type.get(), speaker)

def update_inputs(frame, sim_type, speaker: SpeakerType):
    if str(speaker) != sim_type:
        if sim_type == 'Bass Reflex':
            speaker = BassReflex(UNIT)
        elif sim_type == 'Passive Slave':
            speaker = PassiveSlave(UNIT)
        elif sim_type == '6th Order Bandpass':
            speaker = Bandpass6thOrder(UNIT)

    for widget in frame.winfo_children():
        widget.destroy()

    input_fields = {
        'Bass Reflex': ["Cabinet Volume (L):", "Port Radius (cm):", "Port length (cm):"],
        'Passive Slave': ["Cabinet Volume (L):", "Passive Radiator Compliance (mm/N):", "Passive Radiator Mass (g):", "Passive Radiator Resistance (Ohm):"],
        '6th Order Bandpass': ["Cabinet volume front chamber (L):", "Cabinet volume rear chamber (L):", "Port radius front chamber (cm):", "Port length front chamber (cm):", "Port radius rear chamber (cm):", "Port length rear chamber (cm):"]
    }

    labels = input_fields[sim_type]
    entries = {label: ttk.Entry(frame) for label in labels}

    for label, entry in entries.items():
        ttk.Label(frame, text=label).pack()
        entry.pack()

    frame.entries = entries  # Storing entries in the frame for access during submit

    # Pre-fill with default values
    if sim_type == 'Bass Reflex':
        entries[labels[0]].insert(0, str(speaker.cabinet.volume))
        entries[labels[1]].insert(0, str(speaker.port.radius))
        entries[labels[2]].insert(0, str(speaker.port.length))
    elif sim_type == 'Passive Slave':
        entries[labels[0]].insert(0, str(speaker.cabinet.volume))
        entries[labels[1]].insert(0, str(speaker.slave.Cas))
        entries[labels[2]].insert(0, str(speaker.slave.Mas))
        entries[labels[3]].insert(0, str(speaker.slave.Ras))
    elif sim_type == '6th Order Bandpass':
        entries[labels[0]].insert(0, str(speaker.front_cabinet.volume))
        entries[labels[1]].insert(0, str(speaker.back_cabinet.volume))
        entries[labels[2]].insert(0, str(speaker.front_port.radius))
        entries[labels[3]].insert(0, str(speaker.front_port.length))
        entries[labels[4]].insert(0, str(speaker.back_ports.radius))
        entries[labels[5]].insert(0, str(speaker.back_ports.length))


def submit(frame, speaker: SpeakerType, ax, canvas):
    entries = frame.entries
    try:
        values = {label: float(entry.get()) for label, entry in entries.items()}
        print(f"Submitted values: {values}")  # Debugging output

        # Update the basreflex object based on new inputs
        # For demonstration, this just assumes all values are correctly ordered and given
        # You will need to adjust according to your basreflex and cabinet attributes

        # Here you should call your simulation function with updated parameters
        # Example: frequency, total_sound_pressure, driver_pressure_level, port_pressure_level = simulate_bass_reflex(basreflex)
        # For now, just a dummy plot update:
        ax.clear()
        ax.plot([1, 2, 3], [1, 4, 9])  # Example plot, replace with your actual plotting code
        ax.set_title('Updated Plot')
        canvas.draw()

    except ValueError as e:
        print("Error: Please enter valid numeric values.", e)

if __name__ == "__main__":
    setup_gui()

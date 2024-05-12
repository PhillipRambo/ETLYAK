import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lib2 import unit_from_blue_planet_parameters, SpeakerType, BassReflex, Cabinet, Port, PassiveUnit, PassiveSlave, SimulationType, Bandpass6thOrder, simulate_6thorderbandpass

UNIT = unit_from_blue_planet_parameters(impedance=4, xmax=8, fres=45, bl=8.7, Le=1.18, Re=3.5, Qms=4.37, Qes=0.49, Qts=0.44, Vas=6.8, Sd=127e-3, Mms=38.64e-3, Cms=0.3e-3, Rms=2.57)

speaker = BassReflex(UNIT, Cabinet(8), Port(2, 10))

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

    fig, ax, canvas = setup_plot_window(plot_window)

    setup_input_controls(window, ax, canvas)

    window.mainloop()

def setup_plot_window(window):
    fig = Figure(figsize=(8, 6), dpi=100)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    return fig, ax, canvas

def setup_input_controls(window, ax, canvas):
    frame = ttk.Frame(window, padding="3 3 12 12")
    frame.pack(fill=tk.BOTH, expand=True)

    simulation_types = ['Bass Reflex', 'Passive Slave', '6th Order Bandpass']
    selected_type = tk.StringVar(value=simulation_types[0])

    ttk.Label(frame, text="Select Simulation Type:").pack()
    dropdown = ttk.Combobox(frame, textvariable=selected_type, values=simulation_types, state="readonly")
    dropdown.pack()

    params_frame = ttk.Frame(frame)
    params_frame.pack(fill=tk.BOTH, expand=True, pady=20)

    submit_button = ttk.Button(frame, text="Submit", command=lambda: submit(params_frame, ax, canvas))
    submit_button.pack()

    dropdown.bind("<<ComboboxSelected>>", lambda event: update_inputs(params_frame, selected_type.get(), ax, canvas))
    update_inputs(params_frame, selected_type.get(), ax, canvas)

def update_inputs(frame, sim_type, ax, canvas):
    global speaker
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
    sliders = {label: ttk.Scale(frame, from_=1, to=100, orient='horizontal', command=lambda x=None, f=frame, a=ax, c=canvas: submit(f, a, c)) for label in labels}

    for label, slider in sliders.items():
        ttk.Label(frame, text=label).pack()
        slider.pack()

    frame.sliders = sliders  # Storing sliders in the frame for access during submit

    # Pre-set sliders to default values
    if sim_type == 'Bass Reflex':
        sliders[labels[0]].set(speaker.cabinet.volume)
        sliders[labels[1]].set(speaker.port.radius)
        sliders[labels[2]].set(speaker.port.length)
    elif sim_type == 'Passive Slave':
        sliders[labels[0]].set(speaker.cabinet.volume)
        sliders[labels[1]].set(speaker.slave.Cas)
        sliders[labels[2]].set(speaker.slave.Mas)
        sliders[labels[3]].set(speaker.slave.Ras)
    elif sim_type == '6th Order Bandpass':
        sliders[labels[0]].set(speaker.front_cabinet.volume)
        sliders[labels[1]].set(speaker.back_cabinet.volume)
        sliders[labels[2]].set(speaker.front_port.radius)
        sliders[labels[3]].set(speaker.front_port.length)
        sliders[labels[4]].set(speaker.back_ports.radius)
        sliders[labels[5]].set(speaker.back_ports.length)

def submit(frame, ax, canvas):
    global speaker
    sliders = frame.sliders
    try:
        values = {label: float(slider.get()) for label, slider in sliders.items()}
        # print(f"Submitted values: {values}")  # Debugging output
        if str(speaker) == '6th Order Bandpass':
            # Update Values
            speaker.front_cabinet.volume = values["Cabinet volume front chamber (L):"]
            speaker.back_cabinet.volume = values["Cabinet volume rear chamber (L):"]
            speaker.front_port.radius = values["Port radius front chamber (cm):"]
            speaker.front_port.length = values["Port length front chamber (cm):"]
            speaker.back_ports.radius = values["Port radius rear chamber (cm):"]
            speaker.back_ports.length = values["Port length rear chamber (cm):"]

            # Plot
            f, splT, splF, splR = simulate_6thorderbandpass(speaker)
            ax.clear()
            # ax.semilogx(f, splF)
            # ax.semilogx(f, splR)
            ax.semilogx(f, splT)
            ax.legend(['Front Chamber', 'Rear Chamber', 'Total'])
            ax.set_xlim([f[0], f[-1]])
            ax.set_title('6th Order Bandpass Simulation')
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel('Sound Pressure Level (dB)')
            canvas.draw()

    except ValueError as e:
        print("Error: Please enter valid numeric values.", e)

if __name__ == "__main__":
    setup_gui()

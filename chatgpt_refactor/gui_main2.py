import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from lib import Cabinet, BassReflex, speaker_from_blue_planet_parameters
from simulation_logic import simulate_bass_reflex

def setup_gui():
    # Main window setup for input controls
    window = tk.Tk()
    window.title("Loudspeaker Parameter Input")
    window.geometry("600x400")
    window.configure(bg='#f0f0f0')

    # Separate window for plotting
    plot_window = tk.Toplevel(window)
    plot_window.title("Loudspeaker Simulation Plot")
    plot_window.geometry("800x600")

    fig, ax, canvas = setup_plot_window(plot_window)
    cabinet, speaker_unit, basreflex = setup_speaker_components()

    setup_input_controls(window, cabinet, basreflex, ax, canvas)

    window.mainloop()

def setup_plot_window(window):
    fig = Figure(figsize=(8, 6), dpi=100)
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    return fig, ax, canvas

def setup_input_controls(window, cabinet, basreflex, ax, canvas):
    frame = ttk.Frame(window, padding="3 3 12 12")
    frame.pack(fill=tk.BOTH, expand=True)

    simulation_types = ['Bass Reflex', 'Passive Slave', '6th Order Bandpass']
    selected_type = tk.StringVar(value=simulation_types[0])

    ttk.Label(frame, text="Select Simulation Type:").pack()
    dropdown = ttk.Combobox(frame, textvariable=selected_type, values=simulation_types, state="readonly")
    dropdown.pack()

    params_frame = ttk.Frame(frame)
    params_frame.pack(fill=tk.BOTH, expand=True, pady=20)

    submit_button = ttk.Button(frame, text="Submit", command=lambda: submit(params_frame, cabinet, basreflex, ax, canvas))
    submit_button.pack()

    dropdown.bind("<<ComboboxSelected>>", lambda event: update_inputs(params_frame, selected_type.get()))

def update_inputs(frame, sim_type):
    for widget in frame.winfo_children():
        widget.destroy()

    input_fields = {
        'Bass Reflex': ["Volume (L):", "Width (cm):", "Height (cm):", "Depth (cm):", "Port Size (cm^2):", "Port N:", "Port Length (cm):"],
        'Passive Slave': ["Volume (L):", "Width (cm):", "Height (cm):", "Depth (cm):", "Slave Radius (cm):", "Slave Mass (g):", "Slave Compliance (mm/N):"],
        '6th Order Bandpass': ["Front Volume (L):", "Back Volume (L):", "Front Port Radius (cm):", "Front Port Length (cm):", "Back Port Radius (cm):", "Back Port Length (cm):"]
    }

    labels = input_fields[sim_type]
    entries = {label: ttk.Entry(frame) for label in labels}
    for label, entry in entries.items():
        ttk.Label(frame, text=label).pack()
        entry.pack()

    frame.entries = entries  # Storing entries in the frame for access during submit

def submit(frame, cabinet, basreflex, ax, canvas):
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

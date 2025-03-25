# Script to calculate the capacitance required to run klystron for x amount of useconds (maybe 10?)
# To be determined is whether the klystron can handle a 10% voltage drop. If not, the capacitance must be increased.
# The script will also plot the discharge curve of the RLC circuit with a guessed resistance and inductance.
# calculates the number of capacitors required and the price of the capacitors.
# andthe price of the capacitors.

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
# C = Q/V <-> Q = C * V <-> V = Q/C
# I = dQ/dt -> Q = I * t
# V = Q/C -> V = I * t / C
# C = (I * t) / V
def si_format(value, unit):
    # If the value is 0, just return 0 with the unit
    if value == 0:
        return f"0 {unit}"
    #Define the prefixes
    prefixes = {
        -12: "p",  # pico
        -9: "n",   # nano
        -6: "μ",   # micro
        -3: "m",   # milli
        0: "",     # base
        3: "k",    # kilo
        6: "M",    # mega
        9: "G",    # giga
        12: "T"    # tera
    }
    # Find the exponent of the value by taking the log10 of the absolute value of the value (trust me)
    # and rounding down to the nearest multiple of 3 to fit with the SI prefixes, also called engineering notation I think
    exponent = math.floor(math.log10(abs(value)) / 3) * 3
    exponent = max(min(exponent, 12), -12)  # Clamp between -12 and 12 because that's what we defined and noone uses the other ones anyway
    
    # Scale the value by the exponent
    scaled_value = value / (10 ** exponent)

    # Finally, return the scaled value with the prefix and the unit
    prefix = prefixes.get(exponent, "")
    return f"{scaled_value:.0f} {prefix}{unit}"

V_0 = 35000 # Initial voltage, the max for the klystron
I_requested = 7.57 # Current requested for the klystron
t_run = 10e-6 # Time to run the klystron for
V_drop = 4000 # Voltage drop allowed
t_run_values = [2e-6, 3e-6, 4e-6, 5e-6, 10e-6, 15e-6, 20e-6, 1e-3]
V_drop_values = [500, 1000, 2000, 4000, 5000]
#capacitor_value = 1700e-12 # 1.7nF
capacitor_value = 0.033e-6 # 1.7nF
capacitor_price = 1400 # Kr
def calculate_capacitance(I_requested, t_run, V_drop):
    C = (I_requested * t_run) / V_drop
    return C
def calculate_num_capacitors(C, capacitor_value):
    num_capacitors = math.ceil(C / capacitor_value)
    return num_capacitors
def calculate_actual_capacitance(num_capacitors, capacitor_value):
    actual_capacitance = num_capacitors * capacitor_value
    return actual_capacitance
def calculate_price(num_capacitors, capacitor_price):
    price = num_capacitors * capacitor_price
    return price
def plot_discharge_curve(V_0, R, L, C, t_final):
    """
    Circuit looks like this:
    +---R---L---+
    |           |
    +---C-------+

    Where the L is an equivalent inductance we guess of wires and klystron
    And R is the reesistance we naively assume of the klystron and wires
    The capacitance is the bank we're calculating
    """
    def rlc_equations(y, t):
        """ diff equations for this series RLC circuit.
        y[0] is the voltage across the capacitor
        y[1] is the current in the circuit
        """
        v, i = y
        dvdt = -i / C
        didt = (v - R * i) / L
        return [dvdt, didt]
    
    # Initial conditions: v(0) = V_0, i(0) = 0
    y0 = [V_0, 0]
    
    # Time points
    t = np.linspace(0, t_final, 100000)
    
    # Solve the system of differential equations
    solution = odeint(rlc_equations, y0, t)
    voltage = solution[:, 0]
    current = solution[:, 1]
    
    # Create the plot
    fig, ax1 = plt.subplots(figsize=(10, 6))
        # Convert time to microseconds for x-axis
    t_us = t * 1e6
    
    # Convert voltage to kilovolts for left y-axis
    voltage_kv = voltage / 1e3

    # Voltage plot (left y-axis) in blue
    ax1.set_xlabel('Time (µs)')
    ax1.set_ylabel('Voltage (kV)', color='blue')
    ax1.plot(t_us, voltage_kv, 'b-', label='Voltage') # Plot voltage
    ax1.tick_params(axis='y', labelcolor='blue') # Set the color of the left y-axis to blue
    ax1.grid(True) # Add grid
    
    # Current plot (right y-axis) in red
    ax2 = ax1.twinx() # Create a second y-axis that shares the same x-axis
    ax2.set_ylabel('Current (A)', color='red')
    ax2.plot(t_us, current, 'r-', label='Current') # Plot current
    ax2.tick_params(axis='y', labelcolor='red') # Set the color of the right y-axis to red
    
    # Add vertical dashed line at the run time
    run_time = t_final/4 * 1e6  # Convert to microseconds and find the 1/4 point
    ax1.axvline(x=run_time, color='green', linestyle='--', linewidth=2, label='Run Time End')
   
    # Add a legend in top right corner
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    # Add circuit parameters as text annotation
    num_caps = calculate_num_capacitors(C, capacitor_value)
    actual_cap = calculate_actual_capacitance(num_caps, capacitor_value)
    price = calculate_price(num_caps, capacitor_price)
    circuit_info = (
        f"Circuit Parameters:\n"
        f"Capacitance: {si_format(C, 'F')}\n"
        f"Resistance: {si_format(R, 'Ω')}\n"
        f"Inductance: {si_format(L, 'H')}\n"
        f"Initial Voltage: {si_format(V_0, 'V')}\n"
        f"Price: {price} kr\n"
        f"Number of Capacitors: {num_caps}\n"
        f"Klystron fire time: {si_format(t_final/4, 's')}"  # Since t_final is 4*t_run
    )
        # Position the text box in the upper right corner with a light background round box
    ax1.text(0.07, 0.25, circuit_info, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    plt.title('RLC Circuit Discharge')
    plt.tight_layout()
    plt.show()

def main():
    print("Klystron Supply Calculator")
    print("This script calculates the capacitance required to run a klystron for a certain amount of time.")
    print(f"| {'t_run (s)':>9} | {'V_drop (V)':>10} | {'Capacitance (F)':>15} | {'Num. Caps':>9} | {'Actual Cap (F)':>14} | {'Price (kr)':>10} |")
    print("-"*86)

    for t_run in t_run_values:
        for V_drop in V_drop_values:
            C = calculate_capacitance(I_requested, t_run, V_drop)
            num_caps = calculate_num_capacitors(C, capacitor_value)
            actual_cap = calculate_actual_capacitance(num_caps, capacitor_value)
            price = calculate_price(num_caps, capacitor_price)
            print(f"| {si_format(t_run, 'S'):>9} | {si_format(V_drop,'V'):>10} | {si_format(C, 'F'):>15} | {num_caps:>9} | {si_format(actual_cap, 'F'):>14} | {price:>10} |")

    # Now that we understand a bit the implications of our parameters, we will plot the discharge of a reasonable configuration
    t_run = 10e-6
    V_drop = 2000 # 10% voltage drop
    C = calculate_capacitance(I_requested, t_run, V_drop)
    R = V_0 / I_requested
    L = 100e-6
    t_final = 4 * t_run
    plot_discharge_curve(V_0, R, L, C, t_final)

if __name__ == "__main__":
    main()

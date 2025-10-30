import pyvisa
import time

# Initialize PyVISA resource manager
rm = pyvisa.ResourceManager()

# Open the instrument connection. Replace with your instrument's VISA resource name
psu = rm.open_resource('USB0::0x2A8D::0x1102::MY61002709::INSTR')
psu_txl = rm.open_resource('USB0::0x2A8D::0x1102::MY59004303::INSTR') #probe

# Optional: Identify the instrument
#print(psu.query("*IDN?"))

# Function to ramp voltage from 0 to 3V over 1 second
def ramp_voltage(psu, output, target_voltage=3, ramp_time=1):
    steps = 100  # Number of ramp steps (more steps make it smoother)
    step_voltage = target_voltage / steps
    step_time = ramp_time / steps
    
    # Turn on the output channel
    psu.write(f"OUTP ON, (@{output})")
    
    # Ramp the voltage
    for i in range(steps + 1):
        voltage = i * step_voltage
        psu.write(f"VOLT {voltage}, (@{output})")
        time.sleep(step_time)
    
    # Set the voltage to the final target
    psu.write(f"VOLT {target_voltage}, (@{output})")

# Function to set current limits for source and sink currents
def set_current_limits(psu, output, source_current_limit, sink_current_limit):
    # Set maximum source current
    psu.write(f"CURR {source_current_limit}, (@{output})")
    
    # Set maximum sink current
    psu.write(f"CURR:NEG {sink_current_limit}, (@{output})")

# Function to print PSU status (voltage, current, display text)
def print_psu_status(psu, output):
    voltage = psu.query(f"MEAS:VOLT? (@{output})")
    current = psu.query(f"MEAS:CURR? (@{output})")
    try:
        display_text = psu.query("DISP:TEXT?")  # Query the display text, if supported
    except:
        display_text = "Display text not supported by this PSU."
    
    print(f"Output {output} - Voltage: {voltage.strip()} V, Current: {current.strip()} A")
    print(f"Display Text: {display_text.strip()}")

# Select the outputs for VDD_TDIG, VDD_AFE, and VDD_ADC
output_VDD_DIG_3V = 1  # Assuming channel 1 is VDD_TDIG
output_VDD_6V = 2   # Assuming channel 2 is VDD_AFE
output_VDD_3P3V = 3   # Assuming channel 3 is VDD_ADC

psu.write(f"OUTP OFF, (@{output_VDD_DIG_3V})")
psu.write(f"OUTP OFF, (@{output_VDD_6V})")
psu.write(f"OUTP OFF, (@{output_VDD_3P3V})")

output_VDD_TDIG_3V = 1  # Assuming channel 1 is VDD_TDIG
output_VDD_AFE_3V = 2   # Assuming channel 2 is VDD_AFE
#output_VDD_ADC_3V = 3   # Assuming channel 3 is VDD_ADC

psu_txl.write(f"OUTP OFF, (@{output_VDD_TDIG_3V})")
psu_txl.write(f"OUTP OFF, (@{output_VDD_AFE_3V})")
#psu_txl.write(f"OUTP OFF, (@{output_VDD_ADC_3V})")

# Set the current limit for each output channel to 100e-6 A (100 μA)
current_limit_ch1 = 2e-3 # 2mA minimum for ch1?
current_limit_ch2 = 3e-3 # 1mA minimum for other channels?
current_limit_ch3 = 1e-3 # 1mA minimum for other channels?
set_current_limits(psu, output_VDD_DIG_3V, current_limit_ch1, current_limit_ch1)
set_current_limits(psu, output_VDD_6V, current_limit_ch2, current_limit_ch2)
set_current_limits(psu, output_VDD_3P3V, current_limit_ch3, current_limit_ch3)

# Set the current limit for each output channel to 100e-6 A (100 μA)
current_limit_ch1 = 2e-3 # 2mA minimum for ch1?
current_limit_ch2 = 1e-3 # 1mA minimum for other channels?
current_limit_ch3 = 1e-3 # 1mA minimum for other channels?
set_current_limits(psu_txl, output_VDD_TDIG_3V, current_limit_ch1, current_limit_ch1)
set_current_limits(psu_txl, output_VDD_AFE_3V, current_limit_ch2, current_limit_ch2)
#set_current_limits(psu_txl, output_VDD_ADC_3V, current_limit_ch3, current_limit_ch3)

# Set all outputs to 3V and ramp them up in parallel
psu.write(f"OUTP ON, (@{output_VDD_DIG_3V})")
psu.write(f"OUTP ON, (@{output_VDD_6V})")
psu.write(f"OUTP ON, (@{output_VDD_3P3V})")

ramp_time=5
ramp_voltage(psu, output_VDD_6V, 6, ramp_time)   # Ramp VDD_AFE to 3V
ramp_voltage(psu, output_VDD_3P3V, 3.3, ramp_time)   # Ramp VDD_ADC to 3V
ramp_voltage(psu, output_VDD_DIG_3V, 4, ramp_time)  # Ramp VDD_TDIG to 3V

# Set all outputs to 3V and ramp them up in parallel
psu_txl.write(f"OUTP ON, (@{output_VDD_TDIG_3V})")
psu_txl.write(f"OUTP ON, (@{output_VDD_AFE_3V})")
#psu_txl.write(f"OUTP ON, (@{output_VDD_ADC_3V})")


ramp_time=5
ramp_voltage(psu_txl, output_VDD_TDIG_3V, 4, ramp_time)  # Ramp VDD_TDIG to 3V
ramp_voltage(psu_txl, output_VDD_AFE_3V, 3, ramp_time)   # Ramp VDD_AFE to 3V
#ramp_voltage(psu_txl, output_VDD_ADC_3V, 3, ramp_time)   # Ramp VDD_ADC to 3V


time.sleep(5)

# comment out to remove pyvisa errors
# Print PSU status for each channel
#print_psu_status(psu, output_VDD_DIG_3V)
#print_psu_status(psu, output_VDD_6V)
#print_psu_status(psu, output_VDD_3P3V)


#psu.write(f"OUTP {output_VDD_DIG_3V}, OFF")
#psu.write(f"OUTP {output_VDD_6V}, OFF")
#psu.write(f"OUTP {output_VDD_3P3V}, OFF")

# Close the connection to the instrument
psu.close()
psu_txl.close()
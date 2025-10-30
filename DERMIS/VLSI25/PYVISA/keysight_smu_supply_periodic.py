#!python3
# *********************************************************
# Meas
# USAGE: python.exe keysight_rigol_noise_psd_dermis_dt.py <identifier of this test> <# of runs for averaging> 
#
# Hardware setup:
#   - DERMIS testboard
#   - RIGOLDG3061A - connected to CLK
#   - Keysight MXR058A - oscilloscope (frequency extraction on MATLAB) 
#
# Noise PSD (for DT mode)
# Use RigolDG3061A waveform generator frequency output as CLK. Acquire VOUTP_X0 and VOUTN_X0 with zero differential capacitive value
# Inspired by existing Keysight Infiniium code
# Author: Mark Alea mark.alea@kuleuven.be (ESAT-MICAS, KU LEUVEN)
# Date: 20/08/2024
# *********************************************************
import pyvisa
import string
import struct
import sys, getopt
import scipy
import csv
import numpy as np
import math
import time
from datetime import datetime
import matplotlib.pyplot as plt

# Global variables (booleans: 0 = False, 1 = True).
# ---------------------------------------------------------
debug = 1
# Slow Settling time DC with light (booleans: 0 = False, 1 = True).
# ---------------------------------------------------------
light = 1
# CHANNELS TO BE MEASURED (DERMIS BOARD)
VINP_REPLICA = 1
VINN_REPLICA = 2
VOUTP_REPLICA_INAMP = 3
VOUTN_REPLICA_INAMP = 4
VOUTP_X0_INAMP = 5
VOUTN_X0_INAMP = 6
# =========================================================
# Initialize:
# =========================================================
def initialize():
    # Clear status.
    do_command("*CLS")  
    print("CLS\n")
    # Get and display the device's *IDN? string.
    idn_string = do_query_string("*IDN?")
    print("Identification string: '%s'" % idn_string)
    
    # Load the default setup.
    #do_command("*RST")
    
# =========================================================
# Capture:
# =========================================================
def capture(channel,fin,vamp,vofs,fsample,periods,f_lo):

    '''
    # Set probe attenuation factor.
    do_command(":CHANnel1:PROBe 1.0")
    qresult = do_query_string(":CHANnel1:PROBe?")
    print("Channel 1 probe attenuation factor: %s" % qresult)
    
    # Use auto-scale to automatically set up oscilloscope.
    print("Autoscale.")
    do_command(":AUToscale")
    
    # Set trigger mode.
    do_command(":TRIGger:MODE EDGE")
    qresult = do_query_string(":TRIGger:MODE?")
    print("Trigger mode: %s" % qresult)
    
    # Set EDGE trigger parameters.
    do_command(":TRIGger:EDGE:SOURce CHANnel1")
    qresult = do_query_string(":TRIGger:EDGE:SOURce?")
    print("Trigger edge source: %s" % qresult)
    
    do_command(":TRIGger:LEVel CHANnel1,330E-3")
    qresult = do_query_string(":TRIGger:LEVel? CHANnel1")
    print("Trigger level, channel 1: %s" % qresult)
    
    do_command(":TRIGger:EDGE:SLOPe POSitive")
    qresult = do_query_string(":TRIGger:EDGE:SLOPe?")
    print("Trigger edge slope: %s" % qresult)
    
    # Save oscilloscope setup.
    setup_bytes = do_query_ieee_block(":SYSTem:SETup?")
    
    f = open("setup.set", "wb")
    f.write(setup_bytes)
    f.close()
    print("Setup bytes saved: %d" % len(setup_bytes))
    
    # Change oscilloscope settings with individual commands:
    
    # Set vertical scale and offset.
    do_command(":CHANnel1:SCALe 0.1")
    qresult = do_query_number(":CHANnel1:SCALe?")
    print("Channel 1 vertical scale: %f" % qresult)
    
    do_command(":CHANnel1:OFFSet 0.0")
    qresult = do_query_number(":CHANnel1:OFFSet?")
    print("Channel 1 offset: %f" % qresult)
    '''
    
    # Set horizontal scale and offset.
    #do_command(":TIMebase:SCALe 200e-6")   
    do_command(f":TIMebase:SCALe {0.5/fin}")   
    
    qresult = do_query_string(":TIMebase:SCALe?")
    print("Timebase scale: %s" % qresult)
    
    do_command(":TIMebase:POSition 0.0")
    qresult = do_query_string(":TIMebase:POSition?")
    print("Timebase position: %s" % qresult)
    
    '''# Set the acquisition mode.
    do_command(":ACQuire:MODE RTIMe")
    qresult = do_query_string(":ACQuire:MODE?")
    print("Acquire mode: %s" % qresult)
    # Or, set up oscilloscope by loading a previously saved setup.
    setup_bytes = ""
    f = open("setup.set", "rb")
    setup_bytes = f.read()
    f.close()
    do_command_ieee_block(":SYSTem:SETup", setup_bytes)
    print("Setup bytes restored: %d" % len(setup_bytes))
    # Set the desired number of waveform points,
    # and capture an acquisition.
    do_command(":ACQuire:POINts 32000")
    do_command(":DIGitize")
    '''
    
    do_command(":RUN")
    do_command(":SING")       
    points=fsample*(1/f_lo)*periods  #fixed no. of points to get 1 second windows
    do_command(f":ACQuire:POINts {points}")
    do_command(":DIGitize")
    
    # Download the screen image.
    # --------------------------------------------------------
    screen_bytes = do_query_ieee_block(":DISPlay:DATA? PNG")
    # Save display data values to file.
    filename = 'output\ID_' + inst_id + '_ch_' + str(channel) + '_fin_' + str(fin) + '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) +  '_image.png'
    f = open(filename, "wb")
    f.write(screen_bytes)
    f.close()
    print("Screen image written to screen_image.png.")
    # Download waveform data.
    # --------------------------------------------------------    
# =========================================================
# Capture (FIXED timebase):
# =========================================================
def capture_fixed():

    # Set horizontal scale and offset.
    #do_command(":TIMebase:SCALe 200e-6")   
    do_command(f":TIMebase:SCALe {0.5/1}")   
    
    qresult = do_query_string(":TIMebase:SCALe?")
    print("Timebase scale: %s" % qresult)
    
    do_command(":TIMebase:POSition 0.0")
    qresult = do_query_string(":TIMebase:POSition?")
    print("Timebase position: %s" % qresult)
    
# =========================================================
# Setup Keysight B2901A SMU:
# =========================================================    
def set_smu(instr,vdd_net):

    if vdd_net == "afe_adc_tdig":
        volt=3
        volt_prot=4
        ilimit=1e-3 #expected ~8uA; for tdig peaks
    elif vdd_net == "afe_adc_tdig_dig_3v":
        volt=3
        volt_prot=4
        ilimit=10e-3 #expected ~8uA; for tdig peaks        
    elif vdd_net == "afe":    
        volt=3 #3.3
        volt_prot=4
        ilimit=50e-6 
    elif vdd_net == "adc":    
        volt=3
        volt_prot=6
        ilimit=500e-6        
    elif vdd_net == "tdig":    
        volt=3
        volt_prot=4
        ilimit=1e-3    #for tdig peaks     
    elif vdd_net == "dig_3v":    
        volt=3
        volt_prot=4
        ilimit=10e-3
    elif vdd_net == "dig_6v":    
        volt=6
        volt_prot=6.05
        ilimit=50e-3
    else:
        volt=0
        volt_prot=0
        ilimit=0    
        

    do_command_other(instr,":SOUR:FUNC:MODE VOLT")
    do_command_other(instr,f":SOUR:VOLT {volt}")
    do_command_other(instr,f":SENS:VOLT:PROT {volt_prot}") # 500 uA ilimit    
    
    do_command_other(instr,":SENS:FUNC ""CURR""")    
    do_command_other(instr,":SENS:CURR:RANG:AUTO ON")
    do_command_other(instr,":SENS:CURR:APER 2") # 2 seconds averaging
    do_command_other(instr,f":SENS:CURR:PROT {ilimit}") # 500 uA ilimit
    do_command_other(instr,":OUTP ON")    

import time

def set_smu_ramp(instr, vdd_net, ramp_time=5):
    # Define the voltage and current limits based on vdd_net
    if vdd_net == "afe_adc_tdig":
        target_volt = 3
        volt_prot = 4
        ilimit = 1e-3  # expected ~8uA; for tdig peaks
    elif vdd_net == "afe_adc_tdig_dig_3v":
        target_volt = 3
        volt_prot = 4
        ilimit = 10e-3  # expected ~8uA; for tdig peaks        
    elif vdd_net == "afe":
        target_volt = 3
        volt_prot = 4
        ilimit = 50e-6 
    elif vdd_net == "adc":
        target_volt = 3
        volt_prot = 4
        ilimit = 50e-6        
    elif vdd_net == "tdig":
        target_volt = 3
        volt_prot = 4
        ilimit = 1e-3  # for tdig peaks     
    elif vdd_net == "dig_3v":
        target_volt = 3
        volt_prot = 4
        ilimit = 10e-3
    elif vdd_net == "dig_6v":
        target_volt = 6
        volt_prot = 6.05
        ilimit = 50e-3
    else:
        target_volt = 0
        volt_prot = 0
        ilimit = 0    

    # Configure the SMU for voltage output mode
    do_command_other(instr, ":SOUR:FUNC:MODE VOLT")
    do_command_other(instr, f":SENS:VOLT:PROT {volt_prot}")  # Voltage protection
    
    # Set current measurement settings
    do_command_other(instr, ':SENS:FUNC "CURR"')
    do_command_other(instr, ":SENS:CURR:RANG:AUTO ON")
    do_command_other(instr, ":SENS:CURR:APER 2")  # 2 seconds averaging
    do_command_other(instr, f":SENS:CURR:PROT {ilimit}")  # Current protection
    
    # Turn on output
    do_command_other(instr, ":OUTP ON")
    
    # Ramping voltage logic
    initial_volt = 0  # Start from 0V
    steps = 10  # Number of steps for smooth ramp
    step_time = ramp_time / steps  # Time per step
    volt_step = (target_volt - initial_volt) / steps  # Voltage increment per step
    
    # Ramp voltage from 0 to target voltage over ramp_time seconds
    for step in range(steps):
        current_volt = initial_volt + step * volt_step
        do_command_other(instr, f":SOUR:VOLT {current_volt}")  # Set voltage
        time.sleep(step_time)  # Wait for the next step

    # Finally set the voltage to the target value to ensure precision
    do_command_other(instr, f":SOUR:VOLT {target_volt}")

# Example usage:
# set_smu(my_instr, "afe_adc_tdig", ramp_time=5)


# =========================================================
# Measure Iq from Keysight B2901A SMU:
# =========================================================    
def meas_smu(instr,vdd_net):   
    #do_command_other(instr,":MEAS:CURR? (@1)") 
    do_command_other(instr, ":MEAS:CURR?")
    iq = instr.read()
    iq_per_taxel = float(iq)/4
 
    
    print(f"Iq: {iq}")
    print(f"Iq per taxel: {iq_per_taxel}")
    
    return iq    
    
    
# =========================================================
# Setup Keysight MXR058A Oscilloscope
# =========================================================
def setup_wav(fin,fsample,periods,f_lo):
    do_command(f":ACQuire:SRATe:ANAlog {fsample}") #10MHz sample rate
    do_command(":RUN")
    do_command(":SING")           
    #points=fsample*(1/fin)*periods
    points=fsample*(1/f_lo)*periods  #fixed no. of points to get 1 second windows
    do_command(f":ACQuire:POINts {points}")
    do_command(":DIGitize")
    
# =========================================================
# Setup Keysight MXR058A Oscilloscope
# =========================================================
def init_siggen(instr):
    do_command_other(instr,"*CLS") 
    do_command_other(instr,"SYST:BEEP:STAT ON") 
    time.sleep(2)    
    do_command_other(instr,"SYST:BEEP") 
    time.sleep(2)    
    do_command_other(instr,"OUTP:LOAD INF") #high-Z
    time.sleep(2) 
    
def setup_siggen(instr, fin, vamp, vofs):
    do_command_other(instr,"FUNC SQU") #sine wave
    time.sleep(2)
    do_command_other(instr,f"FREQ +{fin}")  
    time.sleep(2)    
    #do_command_other(instr,f"SOUR:FREQ {fin}")  
    do_command_other(instr,f"VOLT +{vamp}") #pk2pk
    time.sleep(2)    
    do_command_other(instr,f"VOLT:OFFS {vofs}")  
    time.sleep(2)    
    do_command_other(instr,"OUTP ON")     


# =========================================================
# Capture TSENSE:
# =========================================================
def capture_tsense(channel,fin, vamp, vofs):
    # Make measurements.
    # --------------------------------------------------------
    do_command(f":MEASure:SOURce CHANnel{channel}")
    qresult = do_query_string(":MEASure:SOURce?")
    print("Measure source: %s" % qresult)

    do_command(":MEASure:VAVerage")
    qresult = do_query_string(":MEASure:VAVerage?")
    print("TSENSE - Measured average on channel %d: %s" % (channel, qresult))
    
    do_command(":MEASure:VRMS DISP,AC")
    qresult_rms = do_query_string(":MEASure:VRMS? DISP,AC")
    print("TSENSE - Measured sigma on channel %d: %s" % (channel, qresult_rms))    
    
    # Save TSENSE data values to CSV file.
    filename = 'output\ID_' + inst_id + '_tsense_ave' + '_fin_' + str(fin) + '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) +  '.csv'
    f = open(filename, "w")
    f.write("%s\n" % qresult)
    f.write("%s\n" % qresult_rms)
    #Infiniium.write(":DISK:SAVE:WAV CHAN2,""waveform_data"",CSV,ON")
    f.close()
    print(f"Waveform format BYTE data written to {filename}.")   
    
# =========================================================
# Analyze:
# =========================================================
def analyze(channel,fin,vamp,vofs):
    # Make measurements.
    # --------------------------------------------------------
    do_command(f":MEASure:SOURce CHANnel{channel}")
    qresult = do_query_string(":MEASure:SOURce?")
    print("Measure source: %s" % qresult)
    do_command(":MEASure:FREQuency")
    qresult = do_query_string(":MEASure:FREQuency?")
    print("Measured frequency on channel %d: %s" % (channel, qresult))
    do_command(":MEASure:VAMPlitude")
    qresult = do_query_string(":MEASure:VAMPlitude?")
    print("Measured vertical amplitude on channel %d: %s" % (channel, qresult))
    
    if channel==1:
        # Download the screen image.
        # --------------------------------------------------------
        screen_bytes = do_query_ieee_block(":DISPlay:DATA? PNG")
        # Save display data values to file.
        filename = 'output\ID_' + inst_id + '_ch_' + str(channel) + '_fin_' + str(fin) + '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) +  '_image.png'
        f = open(filename, "wb")
        f.write(screen_bytes)
        f.close()
        print("Screen image written to screen_image.png.")
        
    # Download waveform data.
    # --------------------------------------------------------    
    # Get the waveform type.
    qresult = do_query_string(":WAVeform:TYPE?")
    print("Waveform type: %s" % qresult)
    # Get the number of waveform points.
    qresult = do_query_string(":WAVeform:POINts?")
    print("Waveform points: %s" % qresult)
    # Set the waveform source.
    do_command(f":WAVeform:SOURce CHANnel{channel}")
    qresult = do_query_string(":WAVeform:SOURce?")
    print("Waveform source: %s" % qresult)
    # Choose the format of the data returned:
    do_command(":WAVeform:FORMat BYTE")
    print("Waveform format: %s" % do_query_string(":WAVeform:FORMat?"))
    # Display the waveform settings from preamble:
    wav_form_dict = {
    0 : "ASCii",
    1 : "BYTE",
    2 : "WORD",
    3 : "LONG",
    4 : "LONGLONG",
    }
    acq_type_dict = {
    1 : "RAW",
    2 : "AVERage",
    3 : "VHIStogram",
    4 : "HHIStogram",
    6 : "INTerpolate",
    10 : "PDETect",
    }
    acq_mode_dict = {
    0 : "RTIMe",
    1 : "ETIMe",
    3 : "PDETect",
    }
    coupling_dict = {
    0 : "AC",
    1 : "DC",
    2 : "DCFIFTY",
    3 : "LFREJECT",
    }
    units_dict = {
    0 : "UNKNOWN",
    1 : "VOLT",
    2 : "SECOND",
    3 : "CONSTANT",
    4 : "AMP",
    5 : "DECIBEL",
    }
    '''preamble_string = do_query_string(":WAVeform:PREamble?")
    (
    wav_form, acq_type, wfmpts, avgcnt, x_increment, x_origin,   
    x_reference, y_increment, y_origin, y_reference, coupling,
    x_display_range, x_display_origin, y_display_range,
    y_display_origin, date, time, frame_model, acq_mode,
    completion, x_units, y_units, max_bw_limit, min_bw_limit
    ) = preamble_string.split(",")
    print("Waveform format: %s" % wav_form_dict[int(wav_form)])
    print("Acquire type: %s" % acq_type_dict[int(acq_type)])
    print("Waveform points desired: %s" % wfmpts)
    print("Waveform average count: %s" % avgcnt)
    print("Waveform X increment: %s" % x_increment)
    print("Waveform X origin: %s" % x_origin)
    print("Waveform X reference: %s" % x_reference) # Always 0.
    print("Waveform Y increment: %s" % y_increment)
    print("Waveform Y origin: %s" % y_origin)
    print("Waveform Y reference: %s" % y_reference) # Always 0.
    print("Coupling: %s" % coupling_dict[int(coupling)])
    print("Waveform X display range: %s" % x_display_range)
    print("Waveform X display origin: %s" % x_display_origin)
    print("Waveform Y display range: %s" % y_display_range)
    print("Waveform Y display origin: %s" % y_display_origin)
    print("Date: %s" % date)
    print("Time: %s" % time)
    print("Frame model #: %s" % frame_model)
    print("Acquire mode: %s" % acq_mode_dict[int(acq_mode)])
    print("Completion pct: %s" % completion)
    print("Waveform X units: %s" % units_dict[int(x_units)])
    print("Waveform Y units: %s" % units_dict[int(y_units)])
    print("Max BW limit: %s" % max_bw_limit)
    print("Min BW limit: %s" % min_bw_limit)
    '''
    # Get numeric values for later calculations.
    x_increment = do_query_number(":WAVeform:XINCrement?")
    x_origin = do_query_number(":WAVeform:XORigin?")
    y_increment = do_query_number(":WAVeform:YINCrement?")
    y_origin = do_query_number(":WAVeform:YORigin?")
    # Get the waveform data.
    do_command(":WAVeform:STReaming OFF")
    sData = do_query_ieee_block(":WAVeform:DATA?")
    # Unpack signed byte data.
    values = struct.unpack("%db" % len(sData), sData)
    print("Number of data values: %d" % len(values))
    
    # Save waveform data values to CSV file.
    #f = open("waveform_data.csv", "w")
    if channel==VOUTP_X0_INAMP:
        filename = 'output\ID_' + inst_id + '_wav_' + 'voutp_x0' + '-' + str(run_id) + '.csv'
    elif channel==VOUTN_X0_INAMP:   
        filename = 'output\ID_' + inst_id + '_wav_' + 'voutn_x0' + '-' + str(run_id) + '.csv'    
    f = open(filename, "w")
    #f.write("time, voltage\n")
    
    time_arr = []
    voltage_arr = []
    for i in range(0, len(values) - 1):
        time_val = x_origin + (i * x_increment)
        voltage = (values[i] * y_increment) + y_origin
        f.write("%E, %f\n" % (time_val, voltage))
        time_arr.append(time_val)
        voltage_arr.append(voltage)
    #Infiniium.write(":DISK:SAVE:WAV CHAN2,""waveform_data"",CSV,ON")
    f.close()
    print(f"Waveform format BYTE data written to {filename}.")   
    
    time_arr = np.array(time_arr)
    voltage_arr = np.array(voltage_arr)
    
    return time_arr, voltage_arr

# =========================================================
# Round off actual frequency to nearest frequency bin
# =========================================================  
def find_indices(freq, fin, tolerance):
    # Create an empty list to store the filtered indices
    filtered_indices = []
    
    # Loop through each element in the freq array
    for index, value in enumerate(freq):
        # Check if the value is within the range fin*0.5 to fin*1.5
        if fin * tolerance < value < fin * (1+tolerance):
            # If it is, add the index to the filtered_indices list
            filtered_indices.append(index)
    
    return filtered_indices

    
    
# =========================================================
# Perform PSD using welch -- NOT ACTUALLY USED (WELCH DONE IN MATLAB)
# =========================================================  
def extract_power_spectrum(channel,voltage_arr,fin,fsample,vamp,vofs,f_lo):

    windowTime = 1/f_lo #1 #(1/fin) #lower BW = fin/10   #fin/10; 
    windowSize = fsample * windowTime
    overlap = windowSize / 2

    '''
    if channel==1:
        filename = 'output\ID_' + inst_id + '_wav_' + 'vout' + '_' + 'fin_' + str(fin) +  '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) + '.csv'
    else:    
        filename = 'output\ID_' + inst_id + '_wav_' + 'vin' + '_' + 'fin_' + str(fin) +  '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) + '.csv'    
    #f = open(filename, "w")
    #wav = f.read()
    
    wav=[]
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
          wav.append(float(row[1]))    
    
    wav = np.array(wav)
    print(wav[:10])   
    [freq, Pxx] = scipy.signal.welch(x=wav, fs=fsample,nperseg=windowSize, return_onesided=True,scaling='spectrum');    
    '''
    [freq, Pxx] = scipy.signal.welch(x=voltage_arr, fs=fsample,nperseg=windowSize, return_onesided=True,scaling='density');
    freq = np.array(freq)
    Pxx = np.array(Pxx)
    
    print(freq[:10])
    print(Pxx[:10])

    if channel==1:
        filename = 'output\ID_' + inst_id + '_psd_' + 'vout' + '_' + 'fin_' + str(fin) +  '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) + '.csv'
    else:    
        filename = 'output\ID_' + inst_id + '_psd_' + 'vin' + '_' + 'fin_' + str(fin) +  '_vpp_' + str(vamp) + '_vofs_' + str(vofs) + '-' + str(run_id) + '.csv'    
    #f = open(filename, "w")
    #wav = f.read()
    
    f = open(filename, "w")
    #for i in range(0, len(freq) - 1):
    #    f.write("%f, %f\n" % (freq[i], Pxx[i]))
    np.savetxt(filename, np.column_stack((freq,Pxx)), delimiter=',', fmt='%s')
    f.close()
    print(f"Power spectrum format BYTE data written to {filename}.")  

# =========================================================
# Send a command and check for errors:
# =========================================================
def do_command(command, hide_params=False):
    if hide_params:
        (header, data) = command.split(" ", 1)
        if debug:
            print("\nCmd = '%s'" % header)
    else:
        if debug:
            print("\nCmd = '%s'" % command)
            
    Infiniium.write("%s" % command)
    
    if hide_params:
        check_instrument_errors(header)
    else:
        check_instrument_errors(command)
# =========================================================
# Send a command and binary values and check for errors:
# =========================================================
def do_command_ieee_block(command, values):
    if debug:
        print("Cmb = '%s'" % command)
    Infiniium.write_binary_values("%s " % command, values, datatype='B')
    check_instrument_errors(command)
# =========================================================
# Send a query, check for errors, return string:
# =========================================================
def do_query_string(query):
    if debug:
        print("Qys = '%s'" % query)
    result = Infiniium.query("%s" % query)
    check_instrument_errors(query)
    return result
# =========================================================
# Send a query, check for errors, return floating-point value:
# =========================================================
def do_query_number(query):
    if debug:
        print("Qyn = '%s'" % query)
    results = Infiniium.query("%s" % query)
    check_instrument_errors(query)
    return float(results)
# =========================================================
# Send a query, check for errors, return binary values:
# =========================================================
def do_query_ieee_block(query):
    if debug:
        print("Qyb = '%s'" % query)
    result = Infiniium.query_binary_values("%s" % query, datatype='s', container=bytes)
    check_instrument_errors(query, exit_on_error=False)
    return result
# =========================================================
# Check for instrument errors:
# =========================================================
def check_instrument_errors(command, exit_on_error=True):
    while True:
        error_string = Infiniium.query(":SYSTem:ERRor? STRing")
        if error_string: # If there is an error string value.
            if error_string.find("0,", 0, 2) == -1: # Not "No error".
                print("ERROR: %s, command: '%s'" % (error_string, command))
                if exit_on_error:
                    print("Exited because of error.")
                    sys.exit(1)
            else: # "No error"
                break
        else: # :SYSTem:ERRor? STRing should always return string.
            print("ERROR: :SYSTem:ERRor? STRing returned nothing, command: '%s'"% command)
            print("Exited because of error.")
            sys.exit(1)

# =========================================================
# Check for instrument errors:
# =========================================================
def check_instrument_errors_other(instr, command, exit_on_error=True):
    while True:
        error_string = instr.query(":SYSTem:ERRor? STRing")
        if error_string: # If there is an error string value.
            if error_string.find("0,", 0, 2) == -1: # Not "No error".
                print("ERROR: %s, command: '%s'" % (error_string, command))
                if exit_on_error:
                    print("Exited because of error.")
                    sys.exit(1)
            else: # "No error"
                break
        else: # :SYSTem:ERRor? STRing should always return string.
            print("ERROR: :SYSTem:ERRor? STRing returned nothing, command: '%s'"% command)
            print("Exited because of error.")
            sys.exit(1)
    
    
# =========================================================
# Send a command and check for errors: (FOR OTHER INSTRUMENTS)
# =========================================================
def do_command_other(instr, command, hide_params=False):
    if hide_params:
        (header, data) = command.split(" ", 1)
        if debug:
            print("\nCmd = '%s'" % header)
    else:
        if debug:
            print("\nCmd = '%s'" % command)
            
    instr.write("%s" % command)
    
    #if hide_params:
    #    check_instrument_errors_other(instr,header)
    #else:
    #    check_instrument_errors_other(instr,command)
# =========================================================
# Send a command and binary values and check for errors:
# =========================================================

    
# =========================================================
# Main program:
# =========================================================

inst_id =  '_dt_'

#-----------------------------

#rm = pyvisa.ResourceManager("C:\\Windows\\System32\\visa64.dll")
#Infiniium = rm.open_resource("TCPIP0::141.121.231.13::hislip0::INSTR")
rm = pyvisa.ResourceManager();
Infiniium = rm.open_resource("USB0::0x2A8D::0x9007::MY61190114::INSTR")
Infiniium.timeout = 200000 # for 7.5Msamples
Infiniium.clear()

RigolDG3061A = rm.open_resource("USB0::0x1AB1::0x0588::DG3G114600735::INSTR")
RigolDG3061A.timeout = 20000

KeysightB2061A = rm.open_resource("USB0::0x0957::0xCD18::MY51143471::INSTR")
KeysightB2061A.timeout = 200000
KeysightB2061A.clear()

set_smu(KeysightB2061A,"adc") # turn on SMU
#set_smu_ramp(KeysightB2061A,"afe",5) # turn on SMU

#init_siggen(RigolDG3061A) # turn on siggen

#time.sleep(20) #allow for everything to settle after running

  


# Initialize the oscilloscope, capture data, and analyze.
initialize()
#capture()
fsample=500e3 #better matches white noise simulation
periods=10*1.5 #for overlapping 10 windows (50%)
fin=1e6  #0.001 #1e6  # out of band since we need siggen on to set vofs
f_lo=1


vamp=3
fclk=5e3
vofs=1.5

    
iq_array=[]
time_array=[]

    

mins=60
start_time=time.time()
end_time=start_time+ mins*60


while time.time() < end_time:

    iq=meas_smu(KeysightB2061A,"adc")  # record Iq
    print(iq)
    #iq_per_taxel = float(iq)/192
    iq_per_taxel = float(iq)/4 #imode
    print(iq_per_taxel)

    iq_array.append(iq_per_taxel)
    time_array.append(time.time() - start_time)
    
    time.sleep(1)

timestamp=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")    
iq_array=np.array(iq_array)
time_array=np.array(time_array)
combined_data=np.column_stack((time_array,iq_array))
# Calculate the average
average_iq = np.mean(iq_array)

# Print the average
print(f"The average IQ is: {average_iq}")

np.savetxt(f"output\iq_array_{timestamp}.csv",combined_data,delimiter=",",header="Time,Iq")    

plt.plot(time_array,iq_array,marker='o')
plt.xlabel('Time (s)')
plt.ylabel('Iq (AVDD+LCVDD)')
plt.grid(True)
plt.show()
    
Infiniium.close()
RigolDG3061A.close()
KeysightB2061A.close()
print("End of program.")
sys.exit()
#'USB0::0x1AB1::0x0588::DG3G114600735::INSTR'    
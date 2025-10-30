# intuitive
All the files used during the INTUITIVE project.

# BioCAS22
 - MATLAB files used in BioCAS22 paper: https://ieeexplore.ieee.org/document/9948546
   - texture image to tactile sensor data
   - sensor to spikes/frames   

 - Used Kylberg image texture datasets transformed into multi-sensor data (8x8 channels) as input to a modelled spike and frame-based tactile sensor readout.
 - Spikes are used for classification using a Spiking Neural Network trained with Surrogate Gradient (from: Ali Safa).
 - Frames are classified using a classic vanilla Recurrent Neural Network.

 - Run -> N_PA_gen_apdsm_tau_small_sigma_ALL_Kylberg_sweep.m

-------------------------------------------------------
# VLSI23
 - MATLAB files used in VLSI23 paper: https://ieeexplore.ieee.org/document/10185346
   - automated dataset recording -> dataset_workflow_onedac.m
     - texture image to tactile sensor data
     - generate scaled datastream as input for VDACs (used as ideal input for the INTUITIVE chip)
     - automated chip output spikes recording (AER) and address remapping (to fill 12x16 addresses)

  - SNN training code (Python) -> intuitivespikes_snn_sg.py (RUN ON GLIESE SERVER)
      - Texture 180mVpp
    #CUDA_VISIBLE_DEVICES=3 python intuitivespikes_snn_sg.py --datasetPath /users/micas/malea/Documents/ArduinoRead_2/ArduinoRead/out/N_4b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_180mvpp_250mvofs_onedac/ --epochs 150 --hidden_size 256 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes
    
    -  Vibration
    #CUDA_VISIBLE_DEVICES=0 python intuitivespikes_snn_sg.py --datasetPath /users/micas/malea/Documents/ArduinoRead_2/ArduinoRead/out/N_3b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_vibration/ --epochs 150 --hidden_size 128 --num_classes 6 --batch_size 32 --tmax_ms 400 --Ts 10 --kfolds 6 --cache yes --save_bestmodel yes
  - Plots -> slayerperf_plot_Kylberg_chipclassifier.m
     - Requires:
       - 'chip_spikeclassifier' -> SNN training results
       - 'datasetPath' -> spike dataset 

-------------------------------------------------------
# TBCAS24 (JOURNAL) 
- published at: https://doi.org/10.1109/TBCAS.2024.3387545
- intuitive_classifier
	- contains result of characterization of effect of circuit non-idealities to classification accuracy
- intuitive_classifier_sensor
	- contains result of classification of indentation period to integrated e-skin chip sensors
- INTUITIVE_plots_csv_forjournal
	- contains scripts and csv files to plot characterization results of e-skin chip circuits
- out
	- contains previous datasets in VLSI23 paper
- output
	- contains link to Google drive -> output from Pyvisa-based chip characterization
- Poling_waveforms
	- contains measurements during poling, poling setup, etc
- PYVISA_measurements
	- contains python scripts to run PYVISA-based chip characterization. Also includes MATLAB scripts to plot/process results
- rasterplot
	- contains MATLAB scripts to show/process spikes
- SENSOR_DATASETS
	- contains spike output dataset from chip characterization
- SNN_SVM_PLOTS
	- contains scripts to plot SNN+SVM training results
- spikegen
	- contains MATLAB scripts to generate model spikes used in characterization of effect of circuit non-idealities
- utils
	- various useful MATLAB scripts.
 - intuitivechip_classifier_lcsvsnlcs_teensy_run2 -> contains LCS vs N-LCS comparisons (using chip spikes)
 	-  SENSOR_DATASETS -> LCS and N-LCS spike datasets;   SNN_SVM_PLOTS -> contains code to plot SNN+SVM results

-------------------------------------------------------
# ISCAS25
- published at: https://doi.org/10.1109/ISCAS56072.2025.11043712
- DERMIS/ISCAS25
	- AD7746_DERMIS_ALLCHANNELS: Arduino code used to control the capacitive-to-digital evalboard
	- arduino_yy_YshearForce1_F1.3_TR167.csv: Measurement results for the capacitive sensor
	- plot_sensor.m: code to plot measurement results

   
# VLSI25
- published at: https://doi.org/10.23919/VLSITechnologyandCir65189.2025.11074978
- DERMIS/VLSI25/for_paper/BODEPLOT
	- process_gain_cmgain_freq_sweep_multisample.m -> code to plot gain and cmgain
- DERMIS/VLSI25/for_paper/NOISe
	- noise_psd_welsch_inband_inputcsv_wav_ave_Cs_referred.m -> code to plot the noise PSD 
	- link_to_data_sample4.txt -> google drive link for the results (files too big to be downloaded)
- DERMIS/VLSI25/for_paper/SPIKES
	- reconstruct_sndr_read_fit_fin_1_vpp_3.m -> code to reconstruct the signal PSD from spikes
- DERMIS/VLSI25/out -> contains ADC output spikes
- DERMIS/VLSI25/output -> contains gain and cmgain data
- DERMIS/VLSI25/PYVISA -> code to control the instruments via Pyvisa     

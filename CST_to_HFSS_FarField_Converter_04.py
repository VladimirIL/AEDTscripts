# Far-Field data converter from FFS (CST) to FFD (HFSS)
# version 0.4
# supports single- and multi-frequency FF datasets generated in CST as FFS-files (version 3.1 and above)
# Developed by V. Litun, Dipl.-Ing.
# The script is partially based on the earlier version developed by colleagues for the previous revision of the FFS format
# In the current version, the scientific format with 14 decimal digits is employed to represent frequency points

import os
import sys
import numpy as np


working_Directory =r'D:\TEST'
input_file_name = 'farfield_CST_test.ffs'
output_file_name = 'farfield_HFSS_test.ffd'

FF_data_proc = 1 # 0 - 'no processing', '1' - 'scaling incident power to 1 W', '2' - 'scaling accepted power to 1 W', '3' - 'scaling radiated power to 1 W'
flag_freq_linspace = True # flag - forced linearization of frequency samples

field_data = []
global no_of_freq_points, freq_pts_Hz, freq_pts_Hz_linspace, no_of_phi_samples, sim_pwr, acc_pwr, rad_pwr, no_of_theta_samples;
global theta_start, phi_start, theta_end, phi_end, total_data, E_Theta_re, E_Theta_im, E_Phi_re, E_Phi_im, pt_theta, pt_phi;


def read_cst_file(input_file):
    global no_of_freq_points, freq_pts_Hz, sim_pwr, acc_pwr, rad_pwr, no_of_phi_samples, no_of_theta_samples;
    global theta_start, phi_start, theta_end, phi_end, total_data, E_Theta_re, E_Theta_im, E_Phi_re, E_Phi_im,pt_theta, pt_phi;


    inp_file = open(os.path.join(working_Directory,input_file), 'r')
   
    if not os.path.isfile(os.path.join(working_Directory,input_file)):
        print("File path does not exist. Exiting...")
        sys.exit()

    data = inp_file.readlines()
    #First read the number of frequency points and then the frequency value
    cnt = 0 #initializing counter
    skip = 0 #initializing amount of strings in the data-file to skip
    read_data_flag = False #initializing the flag for reading E-field data
    flag_initial_data = False #initializing the flag for the input data
    flag_field_data = False #initializing the flag for FF data definition
    total_data = {}
    fr_id = -1 #initializing frequency block index

    
   
    for line in data:
        # skip string if reqiuired (works for data blocks)
        if skip>0:
            skip -= 1
            continue

        # Reading the number of frequency points in the data-file
        if "Frequencies" in line:
            no_of_freq_points = int(data[cnt+1])
            freq_pts_Hz = np.zeros(no_of_freq_points)
            sim_pwr = np.zeros(no_of_freq_points)
            acc_pwr = np.zeros(no_of_freq_points)
            rad_pwr = np.zeros(no_of_freq_points)
            freq_pts_Hz  = np.zeros(no_of_freq_points)
            freq_pts_Hz_linspace  = np.zeros(no_of_freq_points)
            no_of_phi_samples = np.zeros(no_of_freq_points, dtype=int)
            no_of_theta_samples = np.zeros(no_of_freq_points, dtype=int)
            phi_start = np.zeros(no_of_freq_points)
            theta_start = np.zeros(no_of_freq_points)
            phi_end = np.zeros(no_of_freq_points)
            theta_end = np.zeros(no_of_freq_points)

        # Reading Frequency points and related Simulated and Radiated power values
        if "Radiated" in line or "Stimulated" in line:
            for fr_pt in range(no_of_freq_points):
                rad_pwr[fr_pt] = data[cnt+1+fr_pt*5]
                acc_pwr[fr_pt] = data[cnt+2+fr_pt*5]
                sim_pwr[fr_pt] = data[cnt+3+fr_pt*5]
                freq_pts_Hz[fr_pt] = data[cnt+4+fr_pt*5]
            skip = (fr_pt+1)*5
            cnt = cnt+skip         

        # Reading number of samples
        if "samples" in line:
            fr_id = fr_id+1
            temp_list = data[cnt+1].split(" ")
            no_of_phi_samples[fr_id] = int(temp_list[0])
            no_of_theta_samples[fr_id] = int(temp_list[1])

            flag_initial_data = True
            

        # Initializing arrays for E-feild datasets    
        if flag_initial_data and not(flag_field_data):
            E_Theta_re = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            E_Theta_im = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            E_Phi_re = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            E_Phi_im = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            pt_theta = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            pt_phi = np.zeros((no_of_freq_points, max(no_of_phi_samples)*max(no_of_theta_samples)))
            flag_field_data = True
             
        #Now get the Phi and Theta Sample Pixelations   
        if "E_Theta" in line and read_data_flag == False:
            # reading theta and phi ranges for a frequency point

            read_data_flag = True

            temp_list = data[cnt+1].split(" ")
            new_list = []
            for entry in temp_list:
                if entry != '':
                    new_list.append(entry)
            phi_start[fr_id] = float(new_list[0])
            theta_start[fr_id] = float(new_list[1])
            #Now construct end points
            temp_list = data[(cnt)+no_of_phi_samples[fr_id]*no_of_theta_samples[fr_id]].split(" ")
            new_list = []
            for entry in temp_list:
                if entry != '':
                    new_list.append(entry)
            phi_end[fr_id] = float(new_list[0])
            theta_end[fr_id] = float(new_list[1])

        if read_data_flag == True:
            for E_str_id in range(no_of_phi_samples[fr_id]*no_of_theta_samples[fr_id]):
                line_E = data[cnt+1+E_str_id]
                rE_data_CST =  np.fromstring(line_E, dtype=float, sep=' ');
                E_Theta_re[fr_id,E_str_id] = rE_data_CST[2]
                E_Theta_im[fr_id,E_str_id] = rE_data_CST[3]
                E_Phi_re[fr_id,E_str_id] = rE_data_CST[4]
                E_Phi_im[fr_id,E_str_id] = rE_data_CST[5]
                pt_phi[fr_id,E_str_id] = rE_data_CST[0]
                pt_theta[fr_id,E_str_id] = rE_data_CST[1]
            read_data_flag = False
            skip = no_of_phi_samples[fr_id]*no_of_theta_samples[fr_id]
            cnt = cnt+skip
            print("Field data is successfully collected for the frequency point # " + str(int(fr_id+1)) + " of " + str(int(no_of_freq_points)) + "    -     " +"{:e}".format(freq_pts_Hz[fr_id])+" Hz")
        
        cnt = cnt+1



def write_hfss_file(output_file):
    global no_of_freq_points, freq_pts_Hz, freq_pts_Hz_linspace, no_of_phi_samples, sim_pwr, acc_pwr, rad_pwr, no_of_theta_samples, FF_data_proc;
    global theta_start, phi_start, theta_end, phi_end, E_Theta_re, E_Theta_im, E_Phi_re, E_Phi_im,pt_theta, pt_phi;
    
    out_file = open(os.path.join(working_Directory,output_file), 'w')

    # Zero - processing field dara if required
    FF_factor = [np.ones(no_of_freq_points), 1/np.sqrt(sim_pwr), 1/np.sqrt(acc_pwr), 1/np.sqrt(rad_pwr)];
  
    freq_pts_Hz_linspace = np.linspace(freq_pts_Hz[0], freq_pts_Hz[-1], no_of_freq_points);

    if flag_freq_linspace:
        freq_val_Hz = freq_pts_Hz_linspace
    else:
        freq_val_Hz = freq_pts_Hz


    #First write start end and number of theta divisions
    out_file.write(str(theta_start[0])+" "+str(theta_end[0])+" "+str(no_of_theta_samples[0]))
    out_file.write("\n")
    out_file.write(str(phi_start[0])+" "+str(phi_end[0])+" "+str(no_of_phi_samples[0]))
    out_file.write("\n")
    out_file.write("Frequencies"+ " "+ str(no_of_freq_points))
    out_file.write("\n")

    for fr_pt in range(no_of_freq_points):
        out_file.write("Frequency"+" "+"{:.14e}".format(freq_val_Hz[fr_pt]))
        out_file.write("\n")
        for th_id in range(no_of_theta_samples[fr_pt]):
            for phi_id in range(no_of_phi_samples[fr_pt]):
                E_str_id = phi_id*no_of_theta_samples[fr_pt] + th_id;
                ETr = E_Theta_re[fr_pt,E_str_id]*FF_factor[FF_data_proc][fr_pt];
                ETi = E_Theta_im[fr_pt,E_str_id]*FF_factor[FF_data_proc][fr_pt];
                EPr = E_Phi_re[fr_pt,E_str_id]*FF_factor[FF_data_proc][fr_pt];
                EPi = E_Phi_im[fr_pt,E_str_id]*FF_factor[FF_data_proc][fr_pt];

                my_data = "{:.15e}".format(ETr) + " " + "{:.15e}".format(ETi) + " " + "{:.15e}".format(EPr) + " " + "{:.15e}".format(EPi)
                out_file.write(str(my_data))
                out_file.write("\n")
    out_file.close()
             
             
             
    
              
read_cst_file(input_file_name)

write_hfss_file(output_file_name)

    

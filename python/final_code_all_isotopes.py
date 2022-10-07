#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import array
import matplotlib.pyplot as plt
import seaborn as sb
from pathlib import Path
import sys
import os
import shutil

#global variables
standard_Dataframe = 0
standard_row = 0
df_sum = 0
inclination = 0


##########  Class for Co57 Isotope  ##########
class Co57_isotope:
    """
    Class that has 4 inner classes for 12x12 and 8x8 SiPM arrangement
    """
    
    # create a 1st Inner class for 12x12 SiPM
    class Co57_AnalyticalMethodSolution_12x12_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 0
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/Co_57_Folder/realistic_geometry_detectors/12x12_SiPM_20000_events/Co57_20000_event_12x12_txt/")

            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function
            try:
                for i in range(18):
                    print(self.sorted_file_list[0])
                    self.main(self.sorted_file_list[0])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                


            
        def main(self, sorted_file):
            self.sorted_file = sorted_file
            #!!! sets number of columns we will keep
            source = self.data_folder / self.sorted_file
            used_columns = list(range(152))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)
            
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 144:], inplace = True, axis = 1)
            
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :144], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 144)
            arr.resize(12, 12)
            df_sum = pd.DataFrame(arr)
           
            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)

            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 110

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 140
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                
                Calculated_Energy = -1.27720E-1  +  1.25898E-1*counts[0] + 1.54061E-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                Calculated_Energy_List.append(Calculated_Energy)
                    
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=2))
                n += -4.2
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=2))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final

    # create a 2nd Inner class for 12x12 SiPM
    class Co57_DataAnalysis_12x12_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            
            Co57_9700_events_btin_110_140kev = [0.6, -2.1, -1.67, 0.04, 0.15, -0.79, -0.18, -2.02, 0.27, -0.06, -1.23, 1.43, -0.15, -1.71, 0.73, -0.91, -0.65, 0.01, 0.13, -0.41, -0.14, -0.45, -1.28, 1.6, -0.15, -1.26, 1.22, 0.87, -0.59, 1.78, 1.23, -1.07, 0.56, 0.59, -0.15, 1.82]
            Co57_19300_events_btin_110_140kev = [0.24, -1.49, -0.17, -0.06, -0.48, 0.19, 0.01, 0.5, 4.24, 0.39, -1.1, 1.31, -0.45, -0.77, 0.71, 0.25, -0.99, 0.71, 1.04, -0.35, -0.04, 0.03, -0.85, 1.14, 0.87, -0.52, 1.06, -0.04, -1.08, 0.9, 0.6, -0.82, -0.02, 0.07, -0.41, 1.82]
            #define means for ax1,ax2,ax3
            mean1 = np.mean(Co57_9700_events_btin_110_140kev)
            mean2 = np.mean(Co57_19300_events_btin_110_140kev)
            

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co57_9700_events_btin_110_140kev)
            variance2 = np.var(Co57_19300_events_btin_110_140kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(9700)
            mean_of_events2 = np.mean(19300)
          


            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co57_9700_events_btin_110_140kev)
            USL1 = max(Co57_9700_events_btin_110_140kev)
            
            LSL2 = min(Co57_19300_events_btin_110_140kev)
            USL2 = max(Co57_19300_events_btin_110_140kev)
            
           
            
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co57_9700_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('Co57 Source (Energy = 122 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Co57_19300_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('Co57 Source (Energy = 122 keV) \n φ Accuracy from upper SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


           

            # Define x and y values
            x_realistic_12x12_SiPM = [mean_of_events1, mean_of_events2]
            y_realistic_12x12_SiPM = [sigma1, sigma2]
            # Plot a simple line chart without any feature
            ax3.plot(x_realistic_12x12_SiPM, y_realistic_12x12_SiPM)
            ax3.set_ylim(0.5,5) 
            ax3.set_title('Co57 Source / 12x12 SiPM \n Relation Sigma - events')
            ax3.set_xlabel('Number of events')
            ax3.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            plt.show()


    # create a 3rd Inner class for 8x8 SiPM
    class Co57_AnalyticalMethodSolution_8x8_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/Co_57_Folder/realistic_geometry_detectors/Co57_8x8_SiPM_10000_events/Co57_8x8_SiPM_10000_events_txt/")

            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function
            try:
                for i in range(18):
                    print(self.sorted_file_list[18])
                    self.main(self.sorted_file_list[18])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                


            
        def main(self, sorted_file):
            self.sorted_file = sorted_file
            #!!! sets number of columns we will keep
            source = self.data_folder / self.sorted_file
            used_columns = list(range(68))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)
            
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 144:], inplace = True, axis = 1)
            
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :64], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 64)
            arr.resize(8,8)
            df_sum = pd.DataFrame(arr)
           
            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)

            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 110

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 140
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                if (counts[0] > 1000):
                    Calculated_Energy = 0.8301 + 0.08791*counts[0] + 5.828e-8*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                    

                else:
                    Calculated_Energy = 14.061 + 0.08333*counts[0] + 3.213e-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Co57 Source (Energy = 122 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('v (Energy = 122 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
        
        def threshold_1(self,layer):
            """
            Function that finds the threshold and return a clean dataframe, while 
            setting all values < threshold = 0 
            """
            self.layer = layer
            
    ##        var = 0
    ##        for i in layer:        
    ##            max_value = layer.iloc[i].max()        
    ##            var += max_value
            
            threshold = 2000
            
            #setting a threshold to layer: Dataframe
            clean_dataframe = layer[layer > threshold] 
            
            clean_dataframe = clean_dataframe.fillna(0)
            return clean_dataframe
            
            """
            # if you want to print the clean_dataframe
            ax = sns.heatmap(new_dataframe, vmin = threshold, annot = True , linewidths=.5 ,square=True )
        
            plt.title("Dataframe after apply threshold for phi=0",fontsize=12)
            plt.show()
            """
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=2))
                n += -4.2
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=2))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final

    # create a 4rth Inner class
    class Co57_DataAnalysis_8x8_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            
            Co57_9750_events_btin_110_140kev = [-0.08, -0.42, -0.57, 0.01, -0.47, 0.46, -0.02, -0.84, 2.32, -0.65, -0.66, 1.06, 1.02, -1.54, -0.57, -1.13, 0.21, 2.56, -1.73, -0.64, -0.05, -0.4, -0.92, 0.2, -0.29, 0.79, -0.01, 0.62, -0.2, -0.55, -0.32, -0.79, 0.17, -0.26, -1.39, -0.15]
            #define means for ax1,ax2,ax3
            mean1 = np.mean(Co57_9750_events_btin_110_140kev)
            

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co57_9750_events_btin_110_140kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(9750)
          


            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co57_9750_events_btin_110_140kev)
            USL1 = max(Co57_9750_events_btin_110_140kev)
                       
            
            num_bins = 80
           
            fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co57_9750_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('Co57 Source (Energy = 122 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            

           

            # Define x and y values
            x_realistic_8x8_SiPM = [mean_of_events1]
            y_realistic_8x8_SiPM = [sigma1]
            # Plot a simple line chart without any feature
            ax2.plot(x_realistic_8x8_SiPM, y_realistic_8x8_SiPM)
            ax2.set_ylim(0.5,5) 
            ax2.set_title('Co57 Source / 8x8 SiPM \n Relation Sigma - events')
            ax2.set_xlabel('Number of events')
            ax2.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            plt.show()



##########  Class for Cs137 Isotope  ##########
class Cs137_isotope():
    """
    Class that has 4 inner classes for 12x12 and 8x8 SiPM arrangement
    """
    
    # create a 1st Inner class for 12x12 SiPM
    class Cs137_AnalyticalMethodSolution_12x12_SiPM:
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/Cs_137_Folder/realistic_geometry_detector/Cs137_12x12_SiPM_30000_events/Cs137_12x12_SiPM_30000_events_txt/")

            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function
            try:
                for i in range(18):
                    print(self.sorted_file_list[18])
                    self.main(self.sorted_file_list[18])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                


            
        def main(self, sorted_file):
            self.sorted_file = sorted_file
            #!!! sets number of columns we will keep
            source = self.data_folder / self.sorted_file
            used_columns = list(range(152))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)
            
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 144:], inplace = True, axis = 1)
            
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :144], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 144)
            arr.resize(12, 12)
            df_sum = pd.DataFrame(arr)
           
            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)

            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 630

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 700
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                
                Calculated_Energy = -1.27720E-1  +  1.25898E-1*counts[0] + 1.54061E-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                Calculated_Energy_List.append(Calculated_Energy)
                    
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Cs137 Source (Energy = 662 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Cs137 Source (Energy = 662 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
        
        def threshold_1(self,layer):
            """
            Function that finds the threshold and return a clean dataframe, while 
            setting all values < threshold = 0 
            """
            self.layer = layer
            
    ##        var = 0
    ##        for i in layer:        
    ##            max_value = layer.iloc[i].max()        
    ##            var += max_value
            
            threshold = 2000
            
            #setting a threshold to layer: Dataframe
            clean_dataframe = layer[layer > threshold] 
            
            clean_dataframe = clean_dataframe.fillna(0)
            return clean_dataframe
            
            """
            # if you want to print the clean_dataframe
            ax = sns.heatmap(new_dataframe, vmin = threshold, annot = True , linewidths=.5 ,square=True )
        
            plt.title("Dataframe after apply threshold for phi=0",fontsize=12)
            plt.show()
            """
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=2))
                n += -4.2
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=2))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final

    # create a 2nd Inner class for 12x12 SiPM
    class Cs137_DataAnalysis_12x12_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            Cs137_5000_events_btin_630_700kev = [-8.66, 3.61, -24.89, -0.99, -6.68, 8.16, 2.13, -0.35, 5.91, 18.09, 2.2, 4.46, 3.3, -14.18, -19.31, -6.57, -12.52, -3.23, -0.08, -22.3, -5.04, 7.35, 5.92, 3.99, 3.94, 1.72, -6.01, -4.59, -13.21, 0.21, 18.49, 0.87, 0.4, 5.57, 16.73, -13.12]
            Cs137_10000_events_btin_630_700kev = [-7.45, -3.4, -6.82, -10.44, -0.59, 15.82, 0.77, 3.2, 3.8, 5.95, 0.89, 2.65, 2.89, -7.14, -11.95, -4.4, -10.99, -3.3, -1.2, -9.87, 2.62, 3.42, 0.44, 10.18, -2.38, -2.22, 2.05, -3.7, -8.25, -1.89, 10.57, -2.97, 5.1, -1.94, 7.23, 1.88]
            Cs137_15000_events_btin_630_700kev = [-4.99, -2.39, 2.84, 4.22, -2.98, 10.46, -2.35, -1.35, 1.85, 0.07, 0.61, 2.99, 5.38, -6.43, -9.8, 3.33, -9.1, -2.71, -4.38, -5.13, 6.27, -2.44, 0.48, 12.18, -4.54, 2.98, 3.18, -0.79, -6.8, -2.34, 7.32, -0.93, 7.63, 0.48, 3.01, 2.59]
            
            
            #define means for ax1,ax2,ax3
            mean1 = np.mean(Cs137_5000_events_btin_630_700kev)
            mean2 = np.mean(Cs137_10000_events_btin_630_700kev)
            mean3 = np.mean(Cs137_15000_events_btin_630_700kev)
            

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Cs137_5000_events_btin_630_700kev)
            variance2 = np.var(Cs137_10000_events_btin_630_700kev)
            variance3 = np.var(Cs137_15000_events_btin_630_700kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(5000)
            mean_of_events2 = np.mean(10000)
            mean_of_events3 = np.mean(15000)
            

            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Cs137_5000_events_btin_630_700kev)
            USL1 = max(Cs137_5000_events_btin_630_700kev)
            
            LSL2 = min(Cs137_10000_events_btin_630_700kev)
            USL2 = max(Cs137_10000_events_btin_630_700kev)
            
            LSL3 = min(Cs137_15000_events_btin_630_700kev)
            USL3 = max(Cs137_15000_events_btin_630_700kev)

           
            
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Cs137_5000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Cs137_10000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Cs137_15000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                         np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 10 deg')
            ax3.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            

            # Define x and y values
            x_realistic_12x12_SiPM = [mean_of_events1, mean_of_events2, mean_of_events3]
            y_realistic_12x12_SiPM = [sigma1, sigma2, sigma3]
    ##        x_realistic_12x12_SiPM = [3500, 3800, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000]
    ##        y_realistic_12x12_SiPM = [35.28, 16.13, 14.09, 11.87, 11.1, 9.844, 11.7, 10.63, 11.03, 11.3]
            # Plot a simple line chart without any feature
            ax4.plot(x_realistic_12x12_SiPM, y_realistic_12x12_SiPM)
            ax4.set_ylim(5,) 
            ax4.set_title('Cs137 Source / 12x12 SiPM \n Relation Sigma - events')
            ax4.set_xlabel('Number of events')
            ax4.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            plt.show()

    # create a 3rd Inner class for 8x8 SiPM
    class Cs137_AnalyticalMethodSolution_8x8_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/Cs_137_Folder/realistic_geometry_detector/Cs137_8x8_SiPM_10000_events/Cs137_8x8_SiPM_10000_events_txt/")

            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function
            try:
                for i in range(18):
                    print(self.sorted_file_list[18])
                    self.main(self.sorted_file_list[18])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                


            
        def main(self, sorted_file):
            self.sorted_file = sorted_file
            #!!! sets number of columns we will keep
            source = self.data_folder / self.sorted_file
            used_columns = list(range(68))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)
            
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 144:], inplace = True, axis = 1)
            
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :64], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 64)
            arr.resize(8,8)
            df_sum = pd.DataFrame(arr)
           
            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)

            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 630

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 700
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                if (counts[0] > 1000):
                    Calculated_Energy = 0.8301 + 0.08791*counts[0] + 5.828e-8*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                    

                else:
                    Calculated_Energy = 14.061 + 0.08333*counts[0] + 3.213e-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Cs137 Source (Energy = 662 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Cs137 Source (Energy = 662 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
        
        def threshold_1(self,layer):
            """
            Function that finds the threshold and return a clean dataframe, while 
            setting all values < threshold = 0 
            """
            self.layer = layer
            
    ##        var = 0
    ##        for i in layer:        
    ##            max_value = layer.iloc[i].max()        
    ##            var += max_value
            
            threshold = 2000
            
            #setting a threshold to layer: Dataframe
            clean_dataframe = layer[layer > threshold] 
            
            clean_dataframe = clean_dataframe.fillna(0)
            return clean_dataframe
            
            """
            # if you want to print the clean_dataframe
            ax = sns.heatmap(new_dataframe, vmin = threshold, annot = True , linewidths=.5 ,square=True )
        
            plt.title("Dataframe after apply threshold for phi=0",fontsize=12)
            plt.show()
            """
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=2))
                n += -4.2
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=2))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final

    # create a 4rth Inner class for 8x8 SiPM
    class Cs137_DataAnalysis_8x8_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            Cs137_5000_events_btin_630_700kev = [8.76, -7.16, -16.67, -1.43, -8.35, -8.46, -21.53, 4.01, -7.41, 8.75, 10.11, -1.57, 11.13, 5.03, 11.98, -10.71, -4.15, 0.31, -9.12, 1.6, -6.78, 1.67, 7.11, 24.52, 6.03, 2.0, 2.56, -3.8, 0.83, -11.77, -19.59, 13.62, 10.5, 5.96, -1.58, 2.54]         
            
            #define means for ax1,ax2,ax3
            mean1 = np.mean(Cs137_5000_events_btin_630_700kev)
            mean2 = np.mean(Cs137_5000_events_btin_630_700kev)
            mean3 = np.mean(Cs137_5000_events_btin_630_700kev)
            

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Cs137_5000_events_btin_630_700kev)
            variance2 = np.var(Cs137_5000_events_btin_630_700kev)
            variance3 = np.var(Cs137_5000_events_btin_630_700kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(5000)
            mean_of_events2 = np.mean(10000)
            mean_of_events3 = np.mean(10500)
            

            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Cs137_5000_events_btin_630_700kev)
            USL1 = max(Cs137_5000_events_btin_630_700kev)
            
            LSL2 = min(Cs137_5000_events_btin_630_700kev)
            USL2 = max(Cs137_5000_events_btin_630_700kev)
            
            LSL3 = min(Cs137_5000_events_btin_630_700kev)
            USL3 = max(Cs137_5000_events_btin_630_700kev)

           
            
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Cs137_5000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Cs137_5000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Cs137_5000_events_btin_630_700kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                         np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 10 deg')
            ax3.set_title('Cs137 Source (Energy = 662 keV) \n φ Accuracy from upper SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            

            # Define x and y values
            x_realistic_8x8_SiPM = [mean_of_events1, mean_of_events2, mean_of_events3]
            y_realistic_8x8_SiPM = [sigma1, sigma2, sigma3]
    ##        x_realistic_12x12_SiPM = [3500, 3800, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000]
    ##        y_realistic_12x12_SiPM = [35.28, 16.13, 14.09, 11.87, 11.1, 9.844, 11.7, 10.63, 11.03, 11.3]
            # Plot a simple line chart without any feature
            ax4.plot(x_realistic_8x8_SiPM, y_realistic_8x8_SiPM)
            ax4.set_ylim(5,) 
            ax4.set_title('Cs137 Source / 8x8 SiPM \n Relation Sigma - events')
            ax4.set_xlabel('Number of events')
            ax4.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            plt.show()

##########  Class for Co60 Isotope  ##########
class Co60_isotope:
    """
    Class that has 4 inner classes for 12x12 and 8x8 SiPM arrangement
    """

    # create a 1st Inner class for 12x12 SiPM
    class Co60_AnalyticalMethodSolution_12x12_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/Co_60_Folder/realistic_geometry_detector/Co60_12x12_SiPM_50000_events/50000_event_12x12_txt/")

            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function
            try:
                for i in range(18):
                    print(self.sorted_file_list[18])
                    self.main(self.sorted_file_list[18])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                


            
        def main(self, sorted_file):
            self.sorted_file = sorted_file
            #!!! sets number of columns we will keep
            source = self.data_folder / self.sorted_file
            used_columns = list(range(152))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)
            
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 144:], inplace = True, axis = 1)
            
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :144], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 144)
            arr.resize(12, 12)
            df_sum = pd.DataFrame(arr)
           
            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)

            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 1150

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 1350
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                
                Calculated_Energy = -1.27720E-1  +  1.25898E-1*counts[0] + 1.54061E-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                Calculated_Energy_List.append(Calculated_Energy)
                    
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
        
        def threshold_1(self,layer):
            """
            Function that finds the threshold and return a clean dataframe, while 
            setting all values < threshold = 0 
            """
            self.layer = layer
            
    ##        var = 0
    ##        for i in layer:        
    ##            max_value = layer.iloc[i].max()        
    ##            var += max_value
            
            threshold = 2000
            
            #setting a threshold to layer: Dataframe
            clean_dataframe = layer[layer > threshold] 
            
            clean_dataframe = clean_dataframe.fillna(0)
            return clean_dataframe
            
            """
            # if you want to print the clean_dataframe
            ax = sns.heatmap(new_dataframe, vmin = threshold, annot = True , linewidths=.5 ,square=True )
        
            plt.title("Dataframe after apply threshold for phi=0",fontsize=12)
            plt.show()
            """
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=2))
                n += -4.2
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=2))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=2))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final


    # create a 2nd Inner class for 12x12 SiPM
    class Co60_DataAnalysis_12x12_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):

            Co60_3500_events_btin_1150_1350kev = [-17.77,-5,-2.59,28.23,-47.84,4.76,25.37,13.79,17.29,26.75,-30.32, -30.31, 5.15,-46.44,13.32,46.69,-30.31,31.08,24.22,-13.38,15.98,45.31,27.22,17.65,-17.77,-10.83,20.16,27.96,17.81,31.04,6.51,-10.28,-5.94,-2.13,2.68,-2.07]
            Co60_7000_events_btin_1150_1350kev = [-1.01, -7.4, -3.04, 13.22, -7.25, 3.58, 11.04, -2.09, 17.51, 20.89, -6.19, 14.71, 0.21, -10.2, -2.59, 10.53, -21.98, 8.7, 13.15, -23.11, 27.28, 26.21, 7.82, 5.41, -15.21, -8.62, 9.23, 13.73, -5.24, -1.96, 4.83, -15.52, 5.22, 23.14, 4.24, -0.08]
            Co60_10500_events_btin_1150_1350kev = [0.09, -5.16, -4.15, 3.9, -8.9, -1.44, -9.08, -2.18, 10.06, 17.63, -9.43, 19.95, 3.41, -13.79, -4.19, 11.88, -22.64, 8.38, 10.11, -21.83, 23.08, 13.61, 8.55, 9.41, -5.28, -7.75, 6.87, 18.12, -8.25, -7.86, 2.63, -13.16, 9.43, 17.49, 0.99, -3.18]
            Co60_14200_events_btin_1150_1350kev = [-1.95, -7.66, -0.08, 3.79, -3.54, -0.68, -8.48, -0.95, 5.0, 15.04, -10.18, 10.83, -7.08, -8.75, -5.67, 5.3, -16.06, 1.6, 2.48, -20.16, 13.15, 12.49, 8.62, 11.27, 1.21, -14.93, 1.92, 20.63, -8.31, -2.91, 3.08, -7.95, 8.91, 18.47, -12.24, -1.99]
            Co60_17700_events_btin_1150_1350kev = [-6.83, -7.8, 2.52, 4.73, 1.8, 3.25, -5.19, 1.24, 8.99, 12.48, -10.46, 6.94, -2.6, -0.34, -8.51, 8.62, -14.05, 4.16, 0.45, -12.76, 9.0, 18.47, 9.84, 10.08, -0.32, -12.39, 1.24, 15.56, -12.48, -2.91, 5.99, -15.17, 6.84, 17.82, -8.57, 0.02]

            #define means for ax1,ax2,ax3
            mean1 = np.mean(Co60_3500_events_btin_1150_1350kev)
            mean2 = np.mean(Co60_7000_events_btin_1150_1350kev)
            mean3 = np.mean(Co60_10500_events_btin_1150_1350kev)
            mean4 = np.mean(Co60_14200_events_btin_1150_1350kev)
            mean5 = np.mean(Co60_17700_events_btin_1150_1350kev)

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co60_3500_events_btin_1150_1350kev)
            variance2 = np.var(Co60_7000_events_btin_1150_1350kev)
            variance3 = np.var(Co60_10500_events_btin_1150_1350kev)
            variance4 = np.var(Co60_14200_events_btin_1150_1350kev)
            variance5 = np.var(Co60_17700_events_btin_1150_1350kev)

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)
            sigma4 = np.sqrt(variance4)
            sigma5 = np.sqrt(variance5)

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(3500)
            mean_of_events2 = np.mean(7000)
            mean_of_events3 = np.mean(10500)
            mean_of_events4 = np.mean(14200)
            mean_of_events5 = np.mean(17700)


            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co60_3500_events_btin_1150_1350kev)
            USL1 = max(Co60_3500_events_btin_1150_1350kev)
            
            LSL2 = min(Co60_7000_events_btin_1150_1350kev)
            USL2 = max(Co60_7000_events_btin_1150_1350kev)
            
            LSL3 = min(Co60_10500_events_btin_1150_1350kev)
            USL3 = max(Co60_10500_events_btin_1150_1350kev)

            LSL4 = min(Co60_14200_events_btin_1150_1350kev)
            USL4 = max(Co60_14200_events_btin_1150_1350kev)

            LSL5 = min(Co60_17700_events_btin_1150_1350kev)
            USL5 = max(Co60_17700_events_btin_1150_1350kev)
            
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co60_3500_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Co60_7000_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Co60_10500_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                         np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 10 deg')
            ax3.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax4 plot
            n, bins_ax4, patches = ax4.hist(Co60_14200_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y4 = ((1 / (np.sqrt(2 * np.pi) * sigma4)) *
                         np.exp(-0.5 * (1 / sigma4 * (bins_ax4 - mean4))**2))
            ax4.plot(bins_ax4, y4, '--',label='best fit', color='red')
            ax4.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax4.set_xlabel('Δφ[deg]')
            ax4.set_ylabel('Events / 10 deg')
            ax4.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax4.annotate("Entries = {}".format(round(mean_of_events4)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Mean = {0:.4g}".format(mean4), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Variance = {0:.4g}".format(variance4), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Sigma = {0:.4g}".format(sigma4), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax5 plot
            n, bins_ax5, patches = ax5.hist(Co60_17700_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y5 = ((1 / (np.sqrt(2 * np.pi) * sigma5)) *
                         np.exp(-0.5 * (1 / sigma5 * (bins_ax5 - mean5))**2))
            ax5.plot(bins_ax5, y5, '--',label='best fit', color='red')
            ax5.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax5.set_xlabel('Δφ[deg]')
            ax5.set_ylabel('Events / 10 deg')
            ax5.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax5.annotate("Entries = {}".format(round(mean_of_events5)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax5.annotate("Mean = {0:.4g}".format(mean5), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax5.annotate("Variance = {0:.4g}".format(variance5), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax5.annotate("Sigma = {0:.4g}".format(sigma5), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


            # Define x and y values
            x_realistic_12x12_SiPM = [mean_of_events1, mean_of_events2, mean_of_events3, mean_of_events4, mean_of_events5]
            y_realistic_12x12_SiPM = [sigma1, sigma2, sigma3, sigma4, sigma5]
    ##        x_realistic_12x12_SiPM = [3500, 3800, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000]
    ##        y_realistic_12x12_SiPM = [35.28, 16.13, 14.09, 11.87, 11.1, 9.844, 11.7, 10.63, 11.03, 11.3]
            # Plot a simple line chart without any feature
            ax6.plot(x_realistic_12x12_SiPM, y_realistic_12x12_SiPM)
            ax6.set_ylim(5,) 
            ax6.set_title('Co60 Source / 12x12 SiPM \n Relation Sigma - events')
            ax6.set_xlabel('Number of events')
            ax6.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            plt.show()

    # create a 3rd Inner class for 8x8 SiPM
    class Co60_AnalyticalMethodSolution_8x8_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])
        
        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/odysseas/diplomatiki/main_code/Co_60_Folder/realistic_geometry_detector/Co60_8x8_SiPM_30000_events/30000_event_8x8_txt/")
            self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  #sort by number due to lambda function

            try:
                for i in range(18):
                    print(self.sorted_file_list[18])
                    self.main(self.sorted_file_list[18])
                    
            except KeyboardInterrupt:
                print("Keyboard Interruption")
            finally:
                print("Final_list_Δφ ",self.final_list_Δφ)
                




        def main(self, sorted_file):
            self.sorted_file = sorted_file
            
            source = self.data_folder / self.sorted_file
            
            #!!! sets number of columns we will keep
            used_columns = list(range(68))
            input_data = np.genfromtxt(source , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
            
            df_raw_data = pd.DataFrame(input_data)   
            self.df_bottom_layer = df_raw_data.copy()
            df_info = df_raw_data.copy()
            
            #Remove last 4 columns that are information
            self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 65:], inplace = True, axis = 1)
            #We keep a DataFrame with total counts for each event in .txt 
            self.sum_bottom_counts = self.df_bottom_layer.sum(axis=1).to_frame()
            
            #keep only significant info / last columns
            df_info.drop(df_info.iloc[:, :64], inplace = True, axis = 1)

            #final dataframe
            df_bottom_layer_specific_region = self.calculateEnergy()   

            #df_sum is 1-D and has Total counts of all events
            global df_sum
            #df_sum = self.df_bottom_layer.sum(axis=0).to_frame()   #computes phi with 10.000 event 
            df_sum = df_bottom_layer_specific_region.sum(axis=0).to_frame()  #computes phi only with 3770 events 
            #convert and resize the df_sum
            arr = df_sum.values.copy()
            arr.resize(1, 64)
            arr.resize(8, 8)
            df_sum = pd.DataFrame(arr)
            #print(df_sum)

            ##########  Calling Functions ########
            
            df_sum = self.threshold_1(df_sum)  #Set a threshold

            self.plot_varycentre_vector(df_sum,self.phi)
            #self.plot_heatmap(df_sum)
            
            self.phi += 10
            print(self.phi)
            self.sorted_file_list.pop(0)
            
        def plot_varycentre_vector(self,dataframe,phi):
            """
            Function that plot a vector from X,Y that are varycentre
            """
            self.phi = float(phi)
            X_varycentre = self.X_varycentre(dataframe)
            Y_varycentre = self.Y_varycentre(dataframe)
            
            #angle_phi = np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi
            self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi)  ,2)
            print("arxiko self.inclination",self.inclination)
            if self.phi == 0:
                mean_error = 0
                Δφ = self.inclination
            elif self.phi < 180:
                if self.inclination >= 0:
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                else:
                    self.inclination += 180
                    Δφ = round((self.phi - abs(self.inclination)),2)
                    mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

            elif self.phi <= 270 :
                self.inclination += 180
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
            else:
                self.inclination += 360
                Δφ = round((self.phi - abs(self.inclination)),2)
                mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)
                
            self.final_list_Δφ.append(Δφ)
            print('------ SiPM INFO ------')
            print("X_var: ",X_varycentre ,"Y_var: " ,Y_varycentre)
            print("Inclination line: " , self.inclination, "and Mean Error: ", mean_error ,"%","and Δφ: ", Δφ ,"[Deg]")
            Vector = np.array([X_varycentre,Y_varycentre])
            ax = plt.axes()
            ax.arrow(0, 0, X_varycentre, Y_varycentre, head_width=0.5, head_length=0.7, fc='lightblue', ec='black')
            
            plt.xlim(-2, 2)
            plt.ylim(-2, 2)
            plt.title('Vector of varycentre coordinates, phi='+ str(self.phi) +'\n Inclination : '+str(self.inclination),fontsize=10)
            print('------------------')

            #plt.savefig('Vector of varycentre coordinates, phi=90', bbox_inches='tight')
            #plt.show()
            
        def calculateEnergy(self):
            """
            Each event produces many counts to SiPM and this function computes the
            initial Energy that g-ray photon had and plot the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
            #detector has 1024 for spectrum
            num_bins = 1024
            
            #User sets Low Boundary Energy and High Boundary Energy
            #LBE = float(input("Enter low boundary Energy [keV]: "))
            LBE = 1150

            #HBE = float(input("Enter high boundary Energy [keV]: "))
            HBE = 1350
            if (HBE - LBE) < 0 :
                print("You entered something wrong!")
                
            
            for index, counts in self.sum_bottom_counts.iterrows():
                if (counts[0] > 1000):
                    Calculated_Energy = 0.8301 + 0.08791*counts[0] + 5.828e-8*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                    

                else:
                    Calculated_Energy = 14.061 + 0.08333*counts[0] + 3.213e-7*counts[0]*counts[0] # from maximum of SiPM counts distribution
                    Calculated_Energy_List.append(Calculated_Energy)
                if (Calculated_Energy) < LBE or (Calculated_Energy > HBE):
                    Rows_Remove_List.append(index)
                    
     
            
            #remove events = counts 
            self.df_bottom_layer_specific_region.drop(Rows_Remove_List,  inplace=True)

            #remove energies that are not in [LBE , HBE]
            df_Calculated_Energy = pd.DataFrame(Calculated_Energy_List)
            df_Calculated_Energy.drop(Rows_Remove_List,  inplace=True)
            
            fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            hist = df_Calculated_Energy.hist(bins = num_bins_ratio,  edgecolor="blue",color='white',ax = ax3)

            #n keeps counts of each bin for for example[ 8. 11.  8. ...  0.  0.  1.]

            #find in which bin belongs the specific energy boundaries
            LBE_bin = int(np.digitize(LBE,bins))
            HBE_bin = int(np.digitize(HBE,bins))

            region_of_interest = HBE_bin - LBE_bin
            #keep counts only in this region of interest [LBE_bin , HBE_bin] = [1000 , 1010]
            del_arr = np.delete(n, np.s_[:LBE_bin], 0)
            sliced_array = del_arr[:region_of_interest + 1]
            #print("Sliced array has length: ",len(sliced_array)," bins")
            print("Between ",LBE, " keV and", HBE," keV, there are ",sum(sliced_array)," events")

            
            ########################
            #setting figure for plot
            
            ax2.set_xlabel('Energy [keV]')
            ax2.set_ylabel('Counts', rotation=90)
            ax2.set_xlim([-1,max(Calculated_Energy_List)+50])
            ax2.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            #plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
            
            
        def plot_heatmap(self,dataframe):
            """
            Functions that prints a heatmap for the specific dataframe
            """    
            #plot a heatmap
            fig, ax = plt.subplots(figsize=(11, 9))
            sb.heatmap(dataframe , vmin = 200)
            plt.show()
        
        
        def threshold_1(self,layer):
            """
            Function that finds the threshold and return a clean dataframe, while 
            setting all values < threshold = 0 
            """
            self.layer = layer
            
    ##        var = 0
    ##        for i in layer:        
    ##            max_value = layer.iloc[i].max()        
    ##            var += max_value
            
            threshold = 2000
            
            #setting a threshold to layer: Dataframe
            clean_dataframe = layer[layer > threshold] 
            
            clean_dataframe = clean_dataframe.fillna(0)
            return clean_dataframe
            
            """
            # if you want to print the clean_dataframe
            ax = sns.heatmap(new_dataframe, vmin = threshold, annot = True , linewidths=.5 ,square=True )
        
            plt.title("Dataframe after apply threshold for phi=0",fontsize=12)
            plt.show()
            """
            
       
            
        def X_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Xi = x coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            df_sum = dataframe.sum().to_frame() #creating 1-d dataframe
            ar = np.array(df_sum)
            ar.resize(1,columns)  #resize  array to (1,8)
            self.all_counts = np.sum(ar)  #find sum number of counts

            print('------ INFO ------')
            print("All counts in dataframe for Co_60 " , round(self.all_counts))
            
            #creating lists for x_coo[] and y_coo[]
            i = dataframe.shape[1]  #i exei tin timi ton sinolikwn columns
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
          
            X_var = sum_var / self.all_counts
            return X_var
            
        
        def Y_varycentre(self,dataframe):
            """
            Function that finds the X varycentre = Sum(Ci*Xi)/Sum(Ci)
            Ci = counter of the specific PMT
            Yi = y coordinate of the PMT ε[6*(-3) , 6*(+3)]
            """
            columns = dataframe.shape[1]
            rows = dataframe.shape[0]
            dataframe_sum = dataframe.sum(axis=1).to_frame() #tora vriskei to athroisma ws pros axona y
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #kanw resize ton array
            all_counts = np.sum(ar) #vrisko ton diaireti mou
            
            
            #dimiourgo tin lista 
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            #i=10
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  # pinakas exei pol/stei me to evros tou x
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195, -3.065, 3.065,
            9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=3.065
            n =-3.065
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(j)
                j+=6.13
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(n)
                n += -6.13
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195, 3.065,
            -3.065, -9.195, -15.325, -21.455]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=3.065
            n = -3.065
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(j)
                j+=6.13
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(n)
                n += -6.13
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final

    # create a 4rth Inner class for 8x8 SiPM
    class Co60_DataAnalysis_8x8_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            Co60_3700_events_btin_1150_1350kev = [11.11, 9.8, 22.84, -17.61, 4.7, 3.93, -21.85, 7.97, -1.65, -20.0, -28.05, -0.07, -2.79, 31.43, 14.84, 0.29, 9.93, 14.19, 16.19, 5.72, 4.76, -2.72, -9.35, -62.34, 13.19, 13.75, -3.24, -25.1, -19.37, 6.37, -30.93, -17.36, -7.16, -14.87, -16.99, 1.4]
            Co60_7500_events_btin_1150_1350kev = [-2.26, 0.88, 10.91, -18.57, -5.4, -0.1, -8.56, 12.08, 1.9, -6.88, -9.69, -5.41, -8.81, 2.41, 4.97, 2.97, 16.44, 24.18, 13.18, -0.71, -10.56, -6.48, -39.3, -19.38, 9.67, 4.3, -4.14, -22.07, -9.86, 7.68, -15.44, 6.44, -5.09, -12.4, -12.98, 8.7]
            Co60_11000_events_btin_1150_1350kev = [-8.37, 7.33, 12.31, -24.14, -7.55, 5.8, -11.38, 14.23, 9.34, 2.67, -13.7, -1.55, 3.52, 4.91, 1.36, -1.91, 12.09, 19.34, 19.32, 1.8, 0.22, 6.08, -36.1, -7.15, 2.3, 3.28, 1.51, -4,8, 9.61, 7.12, -10.13, 8.41, -4.99, -9.77, -7.93, 8.35]

            #define means for ax1,ax2,ax3
            mean1 = np.mean(Co60_3700_events_btin_1150_1350kev)
            mean2 = np.mean(Co60_7500_events_btin_1150_1350kev)
            mean3 = np.mean(Co60_11000_events_btin_1150_1350kev)

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co60_3700_events_btin_1150_1350kev)
            variance2 = np.var(Co60_7500_events_btin_1150_1350kev)
            variance3 = np.var(Co60_11000_events_btin_1150_1350kev)

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(3700)
            mean_of_events2 = np.mean(7500)
            mean_of_events3 = np.mean(11000)

            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co60_3700_events_btin_1150_1350kev)
            USL1 = max(Co60_3700_events_btin_1150_1350kev)
            LSL2 = min(Co60_7500_events_btin_1150_1350kev)
            USL2 = max(Co60_7500_events_btin_1150_1350kev)
            LSL3 = min(Co60_11000_events_btin_1150_1350kev)
            USL3 = max(Co60_11000_events_btin_1150_1350kev)
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co60_3700_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 1.5deg')
            ax1.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Co60_7500_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 1.5deg')
            ax2.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Co60_11000_events_btin_1150_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                         np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-200,200])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 1.5deg')
            ax3.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) \n φ Accuracy from upper SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            # Define x and y values
            x_realistic_8x8_SiPM = [mean_of_events1, mean_of_events2, mean_of_events3]
            y_realistic_8x8_SiPM = [sigma1, sigma2, sigma3]
    ##        x_realistic_12x12_SiPM = [3500, 3800, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000]
    ##        y_realistic_12x12_SiPM = [35.28, 16.13, 14.09, 11.87, 11.1, 9.844, 11.7, 10.63, 11.03, 11.3]
            # Plot a simple line chart without any feature
            ax4.plot(x_realistic_8x8_SiPM, y_realistic_8x8_SiPM)
            ax4.set_ylim(5,) 
            ax4.set_title('Co60 Source / 8x8 SiPM \n Relation Sigma - events')
            ax4.set_xlabel('Number of events')
            ax4.set_ylabel('Sigma')
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.44)
            plt.show()


obj = Co60_isotope()
obj.Co60_DataAnalysis_8x8_SiPM()



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
class AnalyticalMethodSolution:
    
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


class DataAnalysis:
    def __init__(self):
        self.gaussian_plot_distribution()
    
    def gaussian_plot_distribution(self):

        Co60_3500_events_btin_1150_1350kev = [-17.77,-5,-2.59,28.23,-47.84,4.76,25.37,13.79,17.29,26.75,-30.32, -30.31, 5.15,-46.44,13.32,46.69,-30.31,31.08,24.22,-13.38,15.98,45.31,27.22,17.65,-17.77,-10.83,20.16,27.96,17.81,31.04,6.51,-10.28,-5.94,-2.13,2.68,-2.07]
        Co60_7000_events_btin_1150_1350kev = [-1.01, -7.4, -3.04, 13.22, -7.25, 3.58, 11.04, -2.09, 17.51, 20.89, -6.19, 14.71, 0.21, -10.2, -2.59, 10.53, -21.98, 8.7, 13.15, -23.11, 27.28, 26.21, 7.82, 5.41, -15.21, -8.62, 9.23, 13.73, -5.24, -1.96, 4.83, -15.52, 5.22, 23.14, 4.24, -0.08]
        Co60_10500_events_btin_1150_1350kev = [0.09, -5.16, -4.15, 3.9, -8.9, -1.44, -9.08, -2.18, 10.06, 17.63, -9.43, 19.95, 3.41, -13.79, -4.19, 11.88, -22.64, 8.38, 10.11, -21.83, 23.08, 13.61, 8.55, 9.41, -5.28, -7.75, 6.87, 18.12, -8.25, -7.86, 2.63, -13.16, 9.43, 17.49, 0.99, -3.18]
        Co60_14200_events_btin_1150_1350kev = [-1.95, -7.66, -0.08, 3.79, -3.54, -0.68, -8.48, -0.95, 5.0, 15.04, -10.18, 10.83, -7.08, -8.75, -5.67, 5.3, -16.06, 1.6, 2.48, -20.16, 13.15, 12.49, 8.62, 11.27, 1.21, -14.93, 1.92, 20.63, -8.31, -2.91, 3.08, -7.95, 8.91, 18.47, -12.24, -1.99]
        Co60_17700_events_btin_1150_1350kev = [-6.83, -7.8, 2.52, 4.73, 1.8, 3.25, -5.19, 1.24, 8.99, 12.48, -10.46, 6.94, -2.6, -0.34, -8.51, 8.62, -14.05, 4.16, 0.45, -12.76, 9.0, 18.47, 9.84, 10.08, -0.32, -12.39, 1.24, 15.56, -12.48, -2.91, 5.99, -15.17, 6.84, 17.82, -8.57, 0.02]



        Co60_7000_events_btin_550_1350kev = [0.35,1.71,1.9,17.85,-12.38,12.36,12.04,7.96,11.56,22.53,-3.99,17.62,4.6,-18.46,3.79,10.41,-10.24,21.83,4.86,-19.37,16.37,30,5.06,12.04,-6.49,-10,15.13,9.14,15.21,0.28,-7.43,-3.52,-10.21,4.89,-0.29]
        Co60_10000_events_btin_0_1450kev = [1.89, 0.4, -1.42, 16.95, -16.22, 11.4, 9.28, 8.78, 6.03, 18.04, -3.42, 19.87, 3.42, -17.34, 4.94, 8.9, -10.92, 20.32, 4.79, -23.29, 15.54, 30.17, 1.93, 10.97, -5.62, -2.53, 14.55, -7.28, 9.23, 15.6, 2.82, -4.43, -2.34, -7.06, 2.34, 2.82]
        Co60_4500_events_btin_950_1450kev = [-1.21, -1.94, -1.16, 17.34, -16.71, 9.48, 23.5, 10.91, 11.44, 11.48, -14.69, 9.62, 3.34, -20.38, 11.25, 12.95, -11.13, 20.25,8.23, -15.95, 14.13, 31.95, 4.04, 14.24, -9.82, -5.61, 14.33, 10.91, 5.42, 13.68, 9.57, -8.9, -2.51, 0.62, 0.28, 0.62]
        Co60_3800_events_btin_1100_1450kev = [-10.74, -2.71, -2.98, 22.55, -11.76, 6.16, 24.9, 14.86, 13.28, 14.54, -24.44, 12.68, 5.22, -29.39, 15.75, 38.66, -14.11, 27.21, 16.32, -15.63, 12.43, 41.53, 17.94, 15.83, -19.15, -5.55, 16.07, 21.82, 5.14, 8.67, 7.15, -6.47, 3.56, 1.17, 4.46, -4.63]
        Co60_4000_events_btin_1050_1450kev = [-9.41, -1.59, -2.83, 19.25, -18.87, 9.12, 23.58, 14.49, 8.75, 8.43, -16.69, 11.89, 4.69, -26.35, 15.35, 25.69, -11.56, 25.98,12.66, -16.68, 12.91, 38.24, 12.63, 13.48, -15.38, -6.01, 14.98, 13.88, 3.47, 11.09, 10.44, -9.18, 5.2, 2.44, 2.42, 2.48]
        Co60_5000_events_btin_900_1450kev = [-1.69, -2.09, 2.05, 18.3, -13.02, 9.47, 20.5, 9.75, 9.62, 12.65, -11.64, 12.71, 4.48, -16.88, 6.89, 13.71, -11.67, 17.21,7.83, -16.82, 13.18, 31.22, 2.59, 12.91, -9.19, -7.61, 14.66, 9.65, 5.23, 9.7, 5.89, -9.66, -3.3, -2.93, 1.26, 1.3]
        Co60_6000_events_btin_750_1450kev = [0.59, 1.02, 3.48, 16.71, -10.19, 10.01, 16.07, 10.55, 5.48, 15.24, -3.5, 11.07, 5.61, -14.22, 1.66, 9.28, -10.55, 20.89, 4.46, -17.64, 13.17, 29.12, 4.04, 11.37, -7.2, -7.49, 12.47, 7.5, 6.9, 9.65, 6.09, -4.3, -1.25, -4.09, 3.66, 3.27]
        Co60_8000_events_btin_400_1450kev = [2.05, 0.89, 0.56, 16.96, -11.72, 10.4, 11.94, 9.58, 8.22, 16.4, -3.53, 17.49, 3.89, -14.33, 3.69, 7.42, -10.14, 19.55, 5.48, -21.14, 15.92, 31.67, 3.0, 12.3, -6.0, -4.5, 14.49, 7.94, 9.41, 14.07, 4.6, -4.76, -3.77, -7.18, 3.42, 2.28]
        Co60_9000_events_btin_200_1450kev = [1.52, 1.15, -1.31, 16.81, -15.69, 11.48, 9.64, 8.47, 6.01, 17.58, -3.33, 19.37, 3.76, -17.39, 5.1, 8.05, -10.77, 20.08, 4.46, -22.79, 15.56, 29.97, 2.09, 11.34, -6.06, -2.82, 15.23, 7.19, 10.03, 15.25, 2.57, -3.68, -2.29, -6.77, 2.26, 2.5]

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
        


    

#obj = AnalyticalMethodSolution()
analysis = DataAnalysis()










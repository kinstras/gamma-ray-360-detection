#!/usr/bin/env python3
import numpy as np
import pandas as pd
# Implementation of matplotlib spy function
#import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from array import array
#import tensorflow as tf
import matplotlib.pyplot as plt
#from sklearn import model_selection
#from keras.models import Sequential
#from keras.models import Model
#from keras.layers import Dense
#from sklearn.preprocessing import StandardScaler
import seaborn as sb
from pathlib import Path
import sys


#global variables
standard_Dataframe = 0
standard_row = 0
df_sum = 0
inclination = 0
class AnalyticalMethodSolution:
    
    
    def __init__(self):
        data = [] 
        data_sum = []
        
        
        #user sets the path and .txt file
        data_folder = Path("/home/konstantinos/Desktop/odysseas/diplomatiki/main_code/gauss_dist/10000_Co60_v2/")
        name_txt = input("Enter the name of txt file you want to read: ")

        #!! you have to change phi 
        self.phi = 0
        
        file_name = data_folder / name_txt
        
        #!!! sets number of columns we will keep
        used_columns = list(range(148))
        input_data = np.genfromtxt(file_name , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
        
        df_raw_data = pd.DataFrame(input_data)   
        self.df_bottom_layer = df_raw_data.copy()
        df_info = df_raw_data.copy()
        
        #Remove last 4 columns that are information
        self.df_bottom_layer.drop(self.df_bottom_layer.iloc[:, 145:], inplace = True, axis = 1)
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
	
	
    def plot_varycentre_vector(self,dataframe,phi):
        """
        Function that plot a vector from X,Y that are varycentre
        """
        self.phi = float(phi)
        X_varycentre = self.X_varycentre(dataframe)
        Y_varycentre = self.Y_varycentre(dataframe)
        
        #angle_phi = np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi
        self.inclination = round((np.arctan(Y_varycentre/X_varycentre)* 180 / np.pi) ,2)
        
        #print("inclination is :",self.inclination)
        if self.phi == 0:
            mean_error = 0
            Δφ = self.inclination
        else:
            Δφ = round((self.phi - abs(self.inclination)),2)
            mean_error = round((self.phi - abs(self.inclination))*100 / self.phi ,2)

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
        LBE = 1350

        #HBE = float(input("Enter high boundary Energy [keV]: "))
        HBE = 1450
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

        print(max(n))
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
            list1.append(j)
            j+=4.2
        #creating list with negative values of x
        for b in range(0,k):
            list2.append(n)
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
            list1.append(j)
            j+=4.2
        #creating list with negative values of x
        for b in range(0,k):
            list2.append(n)
            n += -4.2
        list1 = list(reversed(list1))
        list_final = list1 + list2
        return list_final




obj = AnalyticalMethodSolution()











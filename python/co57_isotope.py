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
            self.phi = 180
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("your_path_here")

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
           
            ##########  Calling Functions #######

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
            print("initial self.inclination",self.inclination)
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
            initial Energy that g-ray photon had and plots the spectrum
            """
            Calculated_Energy_List = []
            Rows_interest_region_List = []
            Rows_Remove_List = []
            self.df_bottom_layer_specific_region = self.df_bottom_layer.copy()
          
            #detector has 1024 channels for spectrum
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
            
            fig, (ax2 , ax3) = plt.subplots(1,2,figsize=(10,7))

            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            #num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            num_bins_ratio = 100
            
         
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
            ax3.set_title('Co57 Source (Energy = 122 keV) spectrum,\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
            plt.tight_layout()
            plt.show()
            print('------------------')
            return self.df_bottom_layer_specific_region
        
            
       
            
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
            i = dataframe.shape[1]  
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  
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
            dataframe_sum = dataframe.sum(axis=1).to_frame() #find the sum with respect to the y-axis
            #print("Df_sum for y_varycentre : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns) #resize the array
            all_counts = np.sum(ar) #find the divider 
            
            
            #creating the list
            i = dataframe.shape[0] 
            k = int(i/2)
          
            #initialize variables for list1
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list) 
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-23.1, -18.9, -14.7,
            -10.5, -6.3, -2.1, 2.1, 6.3, 10.5, 14.7, 18.9, 23.1] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n =-2.1
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=3))
                j+=4.200
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=3))
                n += -4.200
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [23.1, 18.9, 14.7,
            10.5, 6.3, 2.1, -2.1, -6.3, -10.5, -14.7, -18.9, -23.1]
            """
            
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=2.1
            n = -2.1
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j, ndigits=3))
                j+=4.2
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n, ndigits=3))
                n += -4.2
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final
            
            

    # create a 2nd Inner class for 12x12 SiPM
    class Co57_DataAnalysis_12x12_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()
        
        def gaussian_plot_distribution(self):
            Co57_3500_events_btin_110_140kev = [0.06, -1.18, 4.91, -2.0, -1.7, 0.29, 0.64, 0.05, -1.05, -1.45, -0.05, -0.1, -0.64, -0.4, -1.63, -0.56, 2.35, -0.02, -0.58, -0.83, 0.6, -1.12, -2.15, 0.34, -0.2, 0.66, 1.84, -1.18, 0.11, 0.95, 0.19, -0.01, 1.39, 1.62, 0.34, 1.89]
            Co57_5000_events_btin_110_140kev = [-2.64, -2.34, 0.62, -2.34, 1.99, 2.36, 1.96, -0.68, 0.26, 1.2, 0.66, -0.23, 1.16, 0.11, 2.44, 0.65, -3.41, -0.95, -0.36, 0.58, 1.27, -0.73, -0.39, 0.04, -0.42, -0.83, 2.03, 0.61, -0.33, 0.87, 1.19, 0.71, -0.17, -2.26, 0.06, 1.83]
            Co57_9700_events_btin_110_140kev = [0.6, -2.1, -1.67, 0.04, 0.15, -0.79, -0.18, -2.02, 0.27, -0.06, -1.23, 1.43, -0.15, -1.71, 0.73, -0.91, -0.65, 0.01, 0.13, -0.41, -0.14, -0.45, -1.28, 1.6, -0.15, -1.26, 1.22, 0.87, -0.59, 1.78, 1.23, -1.07, 0.56, 0.59, -0.15, 1.82]
            Co57_15000_events_btin_110_140kev = [-0.47, -2.19, -0.89, -0.8, 0.79, 0.22, 0.51, -1.57, 0.25, 0.37, -0.66, 0.85, 0.29, -1.08, 1.3, -0.38, -1.58, -0.32, -0.04, -0.09, 0.34, -0.55, -0.99, 1.07, -0.25, -1.11, 1.49, 0.77, -0.5, 1.47, 1.22, -0.46, 0.31, -0.38, -0.07, 1.83]
            Co57_19300_events_btin_110_140kev = [0.24, -1.49, -0.17, -0.06, -0.48, 0.19, 0.01, 0.5, 4.24, 0.39, -1.1, 1.31, -0.45, -0.77, 0.71, 0.25, -0.99, 0.71, 1.04, -0.35, -0.04, 0.03, -0.85, 1.14, 0.87, -0.52, 1.06, -0.04, -1.08, 0.9, 0.6, -0.82, -0.02, 0.07, -0.41, 1.82]
           
          
            #define means for ax1,ax2,ax3
            mean1 = np.mean(Co57_3500_events_btin_110_140kev)
            mean2 = np.mean(Co57_5000_events_btin_110_140kev)
            mean3 = np.mean(Co57_9700_events_btin_110_140kev)
            mean4 = np.mean(Co57_15000_events_btin_110_140kev)


            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co57_3500_events_btin_110_140kev)
            variance2 = np.var(Co57_5000_events_btin_110_140kev)
            variance3 = np.var(Co57_9700_events_btin_110_140kev)
            variance4 = np.var(Co57_15000_events_btin_110_140kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)
            sigma4 = np.sqrt(variance4)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(3500)
            mean_of_events2 = np.mean(5000)
            mean_of_events3 = np.mean(10000)
            mean_of_events4 = np.mean(15000)
          


            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co57_3500_events_btin_110_140kev)
            USL1 = max(Co57_3500_events_btin_110_140kev)
            
            LSL2 = min(Co57_5000_events_btin_110_140kev)
            USL2 = max(Co57_5000_events_btin_110_140kev)
            
           
            
            num_bins = 80
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co57_3500_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Co57_5000_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line 
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                         np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Co57_9700_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-15,15])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 10 deg')
            ax3.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax4 plot
            n, bins_ax4, patches = ax4.hist(Co57_15000_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y4 = ((1 / (np.sqrt(2 * np.pi) * sigma4)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax4 - mean4))**2))
            ax4.plot(bins_ax4, y4, '--',label='best fit', color='red')
            ax4.set_xlim([-15,15])
            #ax.set_ylim([0,0.45])
            ax4.set_xlabel('Δφ[deg]')
            ax4.set_ylabel('Events / 10 deg')
            ax4.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax4.annotate("Entries = {}".format(round(mean_of_events4)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Mean = {0:.4g}".format(mean4), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Variance = {0:.4g}".format(variance4), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Sigma = {0:.4g}".format(sigma4), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')


        
            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            #plt.show()
            
            #plotaro se mia grafiki ola ta gauss-fit
            #plt.figure()
            fig, ax = plt.subplots()
            num_bins = 30

            n, bins_ax, patches = plt.hist(Co57_3500_events_btin_110_140kev, bins = num_bins, range=(-10,10),density=True,  edgecolor="blue",color='white')
            n2, bins_ax2, patches2 = plt.hist(Co57_3500_events_btin_110_140kev, bins = 100, range=(-10,10),density=True, alpha = 0)
           
            

            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax2 - mean1))**2))
            
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                 np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax2 - mean3))**2))
            
            y4 = ((1 / (np.sqrt(2 * np.pi) * sigma4)) *
                 np.exp(-0.5 * (1 / sigma4 * (bins_ax2 - mean4))**2))
            
            # Define each line separately
            line1, = plt.plot(bins_ax2, y1, '--',label='photopeak: 3.500 event', color='black')
            line2, = plt.plot(bins_ax2, y2, '--',label='photopeak: 5.000 event', color='red')
            line3, = plt.plot(bins_ax2, y3, '--',label='photopeak: 10.000 event', color='blue')
            line4, = plt.plot(bins_ax2, y4, '--',label='photopeak: 15.000 event', color='purple')

            #filled area
            ax.fill(bins_ax2, y1,facecolor ='black', alpha = 0.2)
            ax.fill(bins_ax2, y2,facecolor ='red', alpha = 0.2)
            ax.fill(bins_ax2, y3,facecolor ='blue', alpha = 0.2)
            ax.fill(bins_ax2, y4,facecolor ='purple', alpha = 0.2)
            ax.set_xlim([-10,10])
            
            plt.legend(loc='upper right')
            #plt.title("Γράφημα resolution της διάταξης \n 8x8 ανιχνευτών πιριτίου ")
            #plt.title("Συγκεντρωτικό γράφημα 'Accuracy Plots' διάταξης \n 12x12 ανιχνευτών πιριτίου του ισοτόπου $Co^{57}$ ")
            plt.xlabel('Δφ[deg]')
            plt.ylabel('Events / 10 deg')
            #plt.title("Γράφημα resolution της διάταξης 12x12 ανιχνευτών πιριτίου ")
            plt.show()


    # create a 3rd Inner class for 8x8 SiPM
    class Co57_AnalyticalMethodSolution_8x8_SiPM:
    
        def extract_integer(self, filename):
            return int(filename.split('.')[0].split('_')[15])

        def __init__(self):
            data = [] 
            data_sum = []
            self.phi = 90
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("your_path_here")
            #self.sorted_file_list = sorted(os.listdir(self.data_folder), key=self.extract_integer)  ####sort by number due to lambda function

            try:
                for i in range(18):
                    print(self.sorted_file_list[9])
                    self.main(self.sorted_file_list[9])
                    
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

            self.plot_varycentre_vector(df_sum,self.phi)
            self.plot_heatmap(df_sum)
            
            
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
            initial Energy that γ-ray photon had and plot the spectrum
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
            i = dataframe.shape[1]  
            k = int(i/2)
            
            final_list = self.create_list_X(k)
            print("X Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list)  
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
            dataframe_sum = dataframe.sum(axis=1).to_frame()
            #print("Df_sum gia y varycentro : " , df_sum)
            ar = np.array(dataframe_sum)
            ar.resize(1,columns)
            all_counts = np.sum(ar) 
            
            
            #creating the list
            i = dataframe.shape[0]  #i exei tin timi ton sinolikwn rows
            k = int(i/2)
            #initialize variables for list1
            
            final_list = self.create_list_Y(k)
            print("Y Coordinates of the PMTs are: ",final_list)
            
            var = np.sum(ar*final_list) 
            sum_var = np.sum(var)
            Y_var = sum_var / self.all_counts
            return Y_var
        
        def create_list_X(self,k):
            """
            Creates the list of X_coordinates = [-21.455, -15.325, -9.195,
            -3.065, 3.065, 9.195, 15.325, 21.455] that each SiPM

            """ 
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=3.065
            n =-3.065
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=3))
                j+=6.13
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=3))
                n += -6.13
            list2 = list(reversed(list2))
            list_final = list2 + list1
            return list_final
                  
        def create_list_Y(self,k):
            """
            Creates the list of Y_coordinates = [21.455, 15.325, 9.195,
            3.065, -3.065, -9.195, -15.325, -21.455]
            """
            list1=[] 
            list2=[]
            #setting the step == distance between two SiPM
            j=3.065
            n = -3.065
            
            #creating list with positive values of x
            for b in range(0,k):
                list1.append(round(j,ndigits=3))
                j+=6.13
            #creating list with negative values of x
            for b in range(0,k):
                list2.append(round(n,ndigits=3))
                n += -6.13
            list1 = list(reversed(list1))
            list_final = list1 + list2
            return list_final
            

    # create a 4rth Inner class
    class Co57_DataAnalysis_8x8_SiPM:
        def __init__(self):
            self.gaussian_plot_distribution()            
        
        def gaussian_plot_distribution(self):
            Co57_3500_events_btin_110_140kev = [-1.9, 3.16, 2.63, 0.86, -0.81, 0.01, 1.2, -1.42, 1.74, 1.05, -4.07, 2.12, -1.31, -0.41, -1.64, 0.59, -0.35, 1.04, 1.3, -1.37, 1.56, 0.26, -1.83, 1.45, -2.07, -1.5, -1.41, -2.68, 0.38, -2.92, 0.62, -0.39, -1.76, 0.2, 0.22, 0.89]
            Co57_5000_events_btin_110_140kev = [-0.24, -0.92, 0.11, 1.29, -1.28, -0.97, 0.11, 0.14, 0.69, -0.56, -1.38, 0.66, -0.5, -0.76, 0.47, 1.06, -2.12, 0.19, -0.48, -1.56, 0.32, 0.37, -1.08, -0.06, 0.26, -1.69, 0.26, -0.91, 0.12, -0.47, -0.18, -0.21, 0.3, -1.4, -0.8, 1.74]
            Co57_9750_events_btin_110_140kev = [-0.08, -0.42, -0.57, 0.01, -0.47, 0.46, -0.02, -0.84, 2.32, -0.65, -0.66, 1.06, 1.02, -1.54, -0.57, -1.13, 0.21, 2.56, -1.73, -0.64, -0.05, -0.4, -0.92, 0.2, -0.29, 0.79, -0.01, 0.62, -0.2, -0.55, -0.32, -0.79, 0.17, -0.26, -1.39, -0.15]
            Co57_15000_events_btin_110_140kev = [-0.14, -0.59, -0.35, 0.45, -0.75, -0.01, 0.02, -0.5, 1.75, -0.61, -0.9, 0.93, 0.48, -1.26, -0.22, -0.38, -0.62, 1.79, -1.28, -0.95, 0.08, -0.14, -0.97, 0.11, -0.1, -0.06, 0.08, 0.08, -0.1, -0.52, -0.27, -0.58, 0.22, -0.64, -1.19, 0.49]
            
          #define means for ax1,ax2,ax3
            mean1 = np.mean(Co57_3500_events_btin_110_140kev)
            mean2 = np.mean(Co57_5000_events_btin_110_140kev)
            mean3 = np.mean(Co57_9750_events_btin_110_140kev)
            mean4 = np.mean(Co57_15000_events_btin_110_140kev)
            

            #define variance for ax1,ax2,ax3
            variance1 = np.var(Co57_3500_events_btin_110_140kev)
            variance2 = np.var(Co57_5000_events_btin_110_140kev)
            variance3 = np.var(Co57_9750_events_btin_110_140kev)
            variance4 = np.var(Co57_15000_events_btin_110_140kev)
            

            #define sigma for ax1,ax2,ax3
            sigma1 = np.sqrt(variance1)
            sigma2 = np.sqrt(variance2)
            sigma3 = np.sqrt(variance3)
            sigma4 = np.sqrt(variance4)
            

            #define mean_of_events for ax1,ax2,ax3
            mean_of_events1 = np.mean(3500)
            mean_of_events2 = np.mean(5000)
            mean_of_events3 = np.mean(10000)
            mean_of_events4 = np.mean(15000)
          


            #define LSL and USL for ax1,ax2,ax3
            LSL1 = min(Co57_3500_events_btin_110_140kev)
            USL1 = max(Co57_3500_events_btin_110_140kev)
                       
            
            num_bins = 50
           
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(12,8))
            #plt.style.use('seaborn-bright')
            
            #seting ax1 plot
            n, bins_ax1, patches = ax1.hist(Co57_3500_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax1 - mean1))**2))
            ax1.plot(bins_ax1, y1, '--',label='best fit', color='red')
            ax1.set_xlim([-20,20])
            #ax.set_ylim([0,0.45])
            ax1.set_xlabel('Δφ[deg]')
            ax1.set_ylabel('Events / 10 deg')
            ax1.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax1.annotate("Entries = {}".format(round(mean_of_events1)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Mean = {0:.4g}".format(mean1), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Variance = {0:.4g}".format(variance1), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax1.annotate("Sigma = {0:.4g}".format(sigma1), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax2 plot
            n, bins_ax2, patches = ax2.hist(Co57_5000_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                 np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            ax2.plot(bins_ax2, y2, '--',label='best fit', color='red')
            ax2.set_xlim([-15,15])
            #ax.set_ylim([0,0.45])
            ax2.set_xlabel('Δφ[deg]')
            ax2.set_ylabel('Events / 10 deg')
            ax2.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax2.annotate("Entries = {}".format(round(mean_of_events2)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Mean = {0:.4g}".format(mean2), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Variance = {0:.4g}".format(variance2), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax2.annotate("Sigma = {0:.4g}".format(sigma2), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax3 plot
            n, bins_ax3, patches = ax3.hist(Co57_9750_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax3 - mean3))**2))
            ax3.plot(bins_ax3, y3, '--',label='best fit', color='red')
            ax3.set_xlim([-15,15])
            #ax.set_ylim([0,0.45])
            ax3.set_xlabel('Δφ[deg]')
            ax3.set_ylabel('Events / 10 deg')
            ax3.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax3.annotate("Entries = {}".format(round(mean_of_events3)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Mean = {0:.4g}".format(mean3), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Variance = {0:.4g}".format(variance3), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax3.annotate("Sigma = {0:.4g}".format(sigma3), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            #seting ax4 plot
            n, bins_ax4, patches = ax4.hist(Co57_15000_events_btin_110_140kev, bins = num_bins, range=(-20,20),density=True,  edgecolor="blue",color='white')
            # add a 'best fit' line
            y4 = ((1 / (np.sqrt(2 * np.pi) * sigma4)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax4 - mean4))**2))
            ax4.plot(bins_ax4, y4, '--',label='best fit', color='red')
            ax4.set_xlim([-15,15])
            #ax.set_ylim([0,0.45])
            ax4.set_xlabel('Δφ[deg]')
            ax4.set_ylabel('Events / 10 deg')
            ax4.set_title('$Co^{57}$ Source (Energy = 122, 136 keV) \n φ Accuracy from bottom SiPM')
            ax4.annotate("Entries = {}".format(round(mean_of_events4)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Mean = {0:.4g}".format(mean4), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Variance = {0:.4g}".format(variance4), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
            ax4.annotate("Sigma = {0:.4g}".format(sigma4), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')

            # set the spacing between subplots
            plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.7)
            #plt.show()
            
            #plotting in a graphic all gauss-dist
            #plt.figure()
            fig, ax = plt.subplots()
            num_bins = 30
            
            n, bins_ax, patches = plt.hist(Co57_3500_events_btin_110_140kev, bins = num_bins, range=(-10,10),density=True,  edgecolor="blue",color='white')
            n2, bins_ax2, patches2 = plt.hist(Co57_3500_events_btin_110_140kev, bins = 100, range=(-10,10),density=True, alpha = 0)
            
            y1 = ((1 / (np.sqrt(2 * np.pi) * sigma1)) *
                 np.exp(-0.5 * (1 / sigma1 * (bins_ax2 - mean1))**2))
            
            y2 = ((1 / (np.sqrt(2 * np.pi) * sigma2)) *
                 np.exp(-0.5 * (1 / sigma2 * (bins_ax2 - mean2))**2))
            
            y3 = ((1 / (np.sqrt(2 * np.pi) * sigma3)) *
                 np.exp(-0.5 * (1 / sigma3 * (bins_ax2 - mean3))**2))
            
            y4 = ((1 / (np.sqrt(2 * np.pi) * sigma4)) *
                 np.exp(-0.5 * (1 / sigma4 * (bins_ax2 - mean4))**2))
            
            # Define each line separately
            line1, = plt.plot(bins_ax2, y1, '--',label='photopeak: 3.500 event', color='black')
            line2, = plt.plot(bins_ax2, y2, '--',label='photopeak: 5.000 event', color='red')
            line3, = plt.plot(bins_ax2, y3, '--',label='photopeak: 10.000 event', color='blue')
            line4, = plt.plot(bins_ax2, y4, '--',label='photopeak: 15.000 event', color='purple')

            #filled area
            ax.fill(bins_ax2, y1,facecolor ='black', alpha = 0.2)
            ax.fill(bins_ax2, y2,facecolor ='red', alpha = 0.2)
            ax.fill(bins_ax2, y3,facecolor ='blue', alpha = 0.2)
            ax.fill(bins_ax2, y4,facecolor ='purple', alpha = 0.2)
            ax.set_xlim([-10,10])
            
            plt.legend(loc='upper right')
            #plt.title("Γράφημα resolution της διάταξης \n 8x8 ανιχνευτών πιριτίου ")
            #plt.title("Συγκεντρωτικό γράφημα 'Accuracy Plots' διάταξης \n 8x8 ανιχνευτών πιριτίου του ισοτόπου $Co^{57}$ ")
            plt.xlabel('Δφ[deg]')
            plt.ylabel('Events / 10 deg')
            #plt.title("Γράφημα resolution της διάταξης 12x12 ανιχνευτών πιριτίου ")
            plt.show()

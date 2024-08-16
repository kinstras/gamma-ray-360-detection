def radar_plot2():
    
    angle = np.arange(0, 360, 10, dtype=float) * np.pi / 180.0
    angle = np.append(angle, angle[0])
    print(len(angle))
    Co57_3500_events_btin_110_140kev = [-1.9, 3.16, 2.63, 0.86, -0.81, 0.01, 1.2, -1.42, 1.74, 1.05, -4.07, 2.12, -1.31, -0.41, -1.64, 0.59, -0.35, 1.04, 1.3, -1.37, 1.56, 0.26, -1.83, 1.45, -2.07, -1.5, -1.41, -2.68, 0.38, -2.92, 0.62, -0.39, -1.76, 0.2, 0.22, 0.89, -1.9]
    Co57_5000_events_btin_110_140kev = [-0.24, -0.92, 0.11, 1.29, -1.28, -0.97, 0.11, 0.14, 0.69, -0.56, -1.38, 0.66, -0.5, -0.76, 0.47, 1.06, -2.12, 0.19, -0.48, -1.56, 0.32, 0.37, -1.08, -0.06, 0.26, -1.69, 0.26, -0.91, 0.12, -0.47, -0.18, -0.21, 0.3, -1.4, -0.8, 1.74, -0.24]
    Co57_9750_events_btin_110_140kev = [-0.08, -0.42, -0.57, 0.01, -0.47, 0.46, -0.02, -0.84, 2.32, -0.65, -0.66, 1.06, 1.02, -1.54, -0.57, -1.13, 0.21, 2.56, -1.73, -0.64, -0.05, -0.4, -0.92, 0.2, -0.29, 0.79, -0.01, 0.62, -0.2, -0.55, -0.32, -0.79, 0.17, -0.26, -1.39, -0.15, -0.08]
    Co57_15000_events_btin_110_140kev = [-0.14, -0.59, -0.35, 0.45, -0.75, -0.01, 0.02, -0.5, 1.75, -0.61, -0.9, 0.93, 0.48, -1.26, -0.22, -0.38, -0.62, 1.79, -1.28, -0.95, 0.08, -0.14, -0.97, 0.11, -0.1, -0.06, 0.08, 0.08, -0.1, -0.52, -0.27, -0.58, 0.22, -0.64, -1.19, 0.49, -0.14]
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    # Define each line separately
    #line1, = plt.plot(angle, Co57_3500_events_btin_110_140kev, c='r',marker="v",ls='-', label='$Co^{57}$')
    #plt.fill_between(angle, Co57_3500_events_btin_110_140kev, step="pre", alpha=0.4)
##    line2, = plt.plot(angle, Co57_5000_events_btin_110_140kev, c='g',marker=(8,2,0),ls='--', label='$Cs^{137}$')
##    line3, = plt.plot(angle, Co57_9750_events_btin_110_140kev, c='m',marker="o",ls='--',fillstyle='none', label='$Co^{60}$')
##    line4, = plt.plot(angle, Co57_15000_events_btin_110_140kev, label='$Co^{5}$')

    ax.set_title("A line plot on a polar axis", va='bottom')
    plt.show()


class plot_many_sources:
    """
    class that plot multiple spectrum source's

    """
    
    class AnalyticalMethodSolution_8x8_SiPM:
    

        def __init__(self):
            pass
            data = [] 
            data_sum = []
            self.phi = 260
            self.final_list_Δφ = []
            #user sets the path and .txt file
            self.data_folder = Path("/home/konstantinos/Desktop/diplomatiki/main_code/multiple_sources_8x8_SiPM/")
            #file = open("path to the file", "r") #r is for read
            self.sorted_file_list = sorted(os.listdir(self.data_folder))  #sort by number due to lambda function
            print("Files: ", self.sorted_file_list)
            try:
                self.main(self.sorted_file_list[0])
                    
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
    
            print(self.phi)
            self.shading_segment_set1()
            
            self.sorted_file_list.pop(0)
            
        def shading_segment_set1(self):
        
            angle = np.arange(0, 360, 10, dtype=float) * np.pi / 180.0
            angle = np.append(angle, angle[0])
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

            ax.fill_between(
            np.linspace(19.57*np.pi/180, 20.43*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            ax.fill_between(
            np.linspace(88.245*np.pi/180, 91.755*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            ax.fill_between(
            np.linspace(257.555*np.pi/180, 262.445*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            
            #ax.set_title("A line plot on a polar axis", va='bottom')
            plt.show()

        
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
            LBE = float(input("Enter low boundary Energy [keV]: "))
            #LBE = 110

            HBE = float(input("Enter high boundary Energy [keV]: "))
            #HBE = 140
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
            
            #fig, (ax2,ax3) = plt.subplots(1,2,figsize=(10,7))
            fig, ax2 = plt.subplots()
            #ax2 holds spectrum between [0-1369] keV
            n, bins, patches = ax2.hist(Calculated_Energy_List, bins = num_bins,  edgecolor="blue",color='white')
            
            #ax3 holds spectrum between [LBE-HBE] keV, number of bins should change 
            #num_bins_ratio = int(round((num_bins*(HBE - LBE)) / max(bins)))
            num_bins_ratio = 512
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
            ax2.set_title('Multiple Source\'s Spectrum (Energy = 122, 136 keV \n & 662 keV & 1173, 1332 keV) ,\n Number of events: {0} '.format(sum(n)))

            ax3.set_xlabel('Energy [keV]')
            ax3.set_ylabel('Counts', rotation=90)
            ax3.set_title('Region of Interest\'s Spectrum (Energy = 662 keV),\n events: {0}  between [{1} keV , {2} keV]'.format(sum(sliced_array),LBE,HBE))
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

    
def shading_segment_set1():
        
            angle = np.arange(0, 360, 10, dtype=float) * np.pi / 180.0
            angle = np.append(angle, angle[0])
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

            ax.fill_between(
            np.linspace(19.57*np.pi/180, 20.43*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            ax.fill_between(
            np.linspace(88.245*np.pi/180, 91.755*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            ax.fill_between(
            np.linspace(257.555*np.pi/180, 262.445*np.pi/180),  # Go from 0 to pi/2
            0,                          # Fill from radius 0
            1,                          # To radius 1
            )
            
            #ax.set_title("A line plot on a polar axis", va='bottom')
            plt.show()

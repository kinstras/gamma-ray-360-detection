from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
from tabulate import tabulate
import os
import shutil
def gauss(x , mean , variance , sigma):
    constant = 1 / (variance*(np.sqrt(2*np.pi)))
    #c = 10
    #print("Constant is :",constant)
    #constant = 1 / (variance*(np.sqrt(2*np.pi())))
    return constant*(np.exp(-(x-mean)**2)) / (2*sigma**2)

    #return (constant*(np.exp(-(x-mean)**2)/(2*sigma**2)))

def plot():
    source_folder_path = '/home/konstantinos/Desktop/odysseas/diplomatiki/data_files_txt/gauss_dist/Eu152_2000/'
    list = sorted(os.listdir(source_folder_path))
    source = source_folder_path + list[0]
    #arr = [-0.16, 1.41, 0.77, -3.34, 0.84, -0.2, 1.81, -1.33, 0.41, 1.26]
    #Co57_data_1000 = [-0.16,1.05, -0.44, -0.45, -1.17, -2.84, 0.24, -0.62, 0.18, 2.52, 0.58, -2.18, -3.19, 0.22, 2.39, 3.75, 3.6, -0.66, -0.17, -0.07, -4.45, 3.27, 2.18, -1.94, 0.24, -0.64, 0.29, -1.15, -0.32, -0.59, -0.74, 0.52, 3.39, 1.1, -1.35, -0.77, 1.94, 1.53, 1.79, -0.29, 0.01, 0.51, 0.58, -1.6, -2.01, 1.61, -1.9 , -0.06, -0.17, 2.78, -2.42, 0.13, 0.02, -1.04, -1.94, 1.48, 1.24, -2.03, 2.26, 0.41, 1.26]
    Eu152_data_1000 = [0.32,-2.53,1.2,-1.35,-1.62,-2.51,-2.17,0.08,2.43,-0.47,-3.13,2.4,-0.88,-1.7,3.63,-0.31,-5.65,3.15,4.53,2.25,0.59,2.23,-1.29,1.75,-3.92,-0.35,-0.34,-7.07,-2.53,3.53,0.23,-3.8,3.55,0.98,1.08,-1.21,-0.71,-1.55,1.64,4.73,1.99,-1.25,-1.02,1.25,-0.71,0.4,1.59,-0.84,-5.28,0.45,3.1,-5.56,0.97,-5.17,2.99,-1.24,-1.31,0.18,3.22]
    #Eu152_data_2000_60 = [-4.77,0.74,-3.57,-0.64,-1.56,-1.82,0.27,-0.83,-2.66,-0.98,1.83,-0.35,-0.87,-2.65,1.56,-1.65,-1.84,1.02,-0.13,0.79,1.56,-2.02,1.69,-0.2,-3.77,-2.73,2.03,-2.03,-0.22,-1.57,-5.25,-2.56,0.49,0.32,0.08,1.91,0.09,4.21,-0.24,0.45,0.28,-0.06,-2.46,-1.41,-2.32,-0.72,-0.74,2.75,-0.06,-1.13,2.47,-1.02,-1.44,2.48,3.69,-0.48,-2.8,-1.41,4.11,-1.31] 
    Eu152_data_2000_50 = [-4.77,0.74,-3.57,-0.64,-1.56,-1.82,0.27,-0.83,-2.66,-0.98,1.83,-0.35,-0.87,-2.65,1.56,-1.65,-1.84,1.02,-0.13,0.79,1.56,-2.02,1.69,-0.2,-3.77,-2.73,2.03,-2.03,-0.22,-1.57,-5.25,-2.56,0.49,0.32,0.08,1.91,0.09,4.21,-0.24,0.45,0.28,-0.06,-2.46,-1.41,-2.32,-0.72,-0.74,2.75,-0.06,-1.13] 
    Eu152_data_2000_40 = [-4.77,0.74,-3.57,-0.64,-1.56,-1.82,0.27,-0.83,-2.66,-0.98,1.83,-0.35,-0.87,-2.65,1.56,-1.65,-1.84,1.02,-0.13,0.79,1.56,-2.02,1.69,-0.2,-3.77,-2.73,2.03,-2.03,-0.22,-1.57,-5.25,-2.56,0.49,0.32,0.08,1.91,0.09,4.21,-0.24,0.45] 
    Co60_data_10000 = [1.53,-6.08,7.15,-6,1.2,4.66,-18.49,5.81,-3.3,7.57,-12.53,3.26,0.19,10.25,16.84,-2.73,14.22,11.25,4.48,9.09,0.76,-1.57,-1.4,14.93,5.95,8.76,8.99,-12.89,-12.19,5.11,-5.48,-6.21,-3.85,1.68,-10.35,0.92]
    #real_data = [-0.16,1.05, -0.44, -0.45, -1.17, -2.84, 0.24, -0.62, 0.18, 2.52, 0.58, -2.18, -3.19, 0.22, 2.39, 3.75]
    #define constants
    mean = np.mean(Co60_data_10000)
    variance = np.var(Co60_data_10000)
    sigma = np.sqrt(variance)
    LSL = min(Co60_data_10000)
    USL = max(Co60_data_10000)
    num_bins = 80
    # calculate the z-transform
    #z1 = ( LSL - mean ) / sigma
    #z2 = ( USL - mean ) / sigma

    #Process capability
    Cp = (USL - LSL)/ (6*variance)
    
##    x = np.arange(z1, z2, 0.001) # range of x in spec
##    x_all = np.arange(-10, 10, 0.001) # entire range of x, both in and out of spec
##
##
##    y = gauss(x , mean, variance , sigma)
##    y_all = gauss(x_all , mean, variance , sigma)
##
    #Print info table
    #info_table = [['delta_phi_upper'],['Entries' , len(real_data)], ['Mean', mean] ,['Variance', variance] ,['Sigma', sigma]]
    #print(tabulate(info_table,headers = 'firstrow', tablefmt = 'fancy_grid'))


    fig, ax = plt.subplots(1,1,figsize=(10,7))
    plt.style.use('seaborn-bright')
    #ax.plot(x_all,y_all, label='Gaussian Distribution', linewidth = 2, color = 'r')
    #ax.legend(title ='delta_phi_upper')
    n, bins, patches = ax.hist(Co60_data_10000, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
    
    
    # add a 'best fit' line
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mean))**2))
    ax.plot(bins, y, '--',label='best fit', color='red')
    ax.set_xlim([-200,200])
    ax.set_ylim([0,0.45])
    ax.set_xlabel('Δφ[deg]')
    ax.set_ylabel('Events / 1.5deg')
    ax.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) φ Accuracy from upper SiPM')
    plt.annotate("Entries = 10.000", xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Mean = {0:.4g}".format(mean), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Variance = {0:.4g}".format(variance), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Sigma = {0:.4g}".format(sigma), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Process Capability = {0:.4g}".format(Cp), xy=(0, 1), xytext=(12, -60), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("NoM = {0:.4g}".format(len(Co60_data_10000)), xy=(0, 1), xytext=(12, -72), va='top',xycoords='axes fraction', textcoords='offset points')

    plt.show()


plot()




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
    return constant*(np.exp(-(x-mean)**2)) / (2*sigma**2)

def plot():
    Co57_1000_events = [-0.16,1.05, -0.44, -0.45, -1.17, -2.84, 0.24, -0.62, 0.18, 2.52, 0.58, -2.18, -3.19, 0.22, 2.39, 3.75, 3.6, -0.66, -0.17, -0.07, -4.45, 3.27, 2.18, -1.94, 0.24, -0.64, 0.29, -1.15, -0.32, -0.59, -0.74, 0.52, 3.39, 1.1, -1.35, -0.77, 1.94, 1.53, 1.79, -0.29, 0.01, 0.51, 0.58, -1.6, -2.01, 1.61, -1.9 , -0.06, -0.17, 2.78, -2.42, 0.13, 0.02, -1.04, -1.94, 1.48, 1.24, -2.03, 2.26, 0.41, 1.26]
    Eu152_1000_events = [0.32,-2.53,1.2,-1.35,-1.62,-2.51,-2.17,0.08,2.43,-0.47,-3.13,2.4,-0.88,-1.7,3.63,-0.31,-5.65,3.15,4.53,2.25,0.59,2.23,-1.29,1.75,-3.92,-0.35,-0.34,-7.07,-2.53,3.53,0.23,-3.8,3.55,0.98,1.08,-1.21,-0.71,-1.55,1.64,4.73,1.99,-1.25,-1.02,1.25,-0.71,0.4,1.59,-0.84,-5.28,0.45,3.1,-5.56,0.97,-5.17,2.99,-1.24,-1.31,0.18,3.22]
    Co60_10000_events = [1.53,-6.08,7.15,-6,1.2,4.66,-18.49,5.81,-3.3,7.57,-12.53,3.26,0.19,10.25,16.84,-2.73,14.22,11.25,4.48,9.09,0.76,-1.57,-1.4,14.93,5.95,8.76,8.99,-12.89,-12.19,5.11,-5.48,-6.21,-3.85,1.68,-10.35,0.92]
    Co60_20000_events = [-4.14,-2.26,0.88,10.91,-18.57,-5.4,-0.1,-8.56,12.08,1.9,-6.88,-9.69,-5.41,-8.81,2.41,4.97,2.97,9.59,24.18,13.18,-0.71,-10.56,-6.48,-39.3,-19.38,9.67,4.3,-4.14,-22.07,-9.86,7.68,-15.44,6.44,-5.09,-12.4,-12.98,8.7]
    Co60_events_between_1150_1350 = [7476,7430,7404,7294,7190,7256,7401,7347,7395,7283,7333]
    #mine measurements to see how sigma changes
    Δφ_Co60_pure_10000 = [1.53,-6.08,7.15,-6,1.2,4.66,-18.49,5.81,-3.3,7.57,-12.53,3.26,0.19,10.25,16.84,-2.73,14.22,11.25,4.48,9.09,0.76,-1.57,-1.4,14.93,5.95,8.76,8.99,-12.89,-12.19,5.11,-5.48,-6.21,-3.85,1.68,-10.35,0.92]
    Δφ_Co_60_pure_20000 = [-2.31,-9.72,-0.36,-11.46,-4.77,3.74,-5.05,4.56,1.85,-2.82,-6.05,-2.23,-3.93,1.7,5.67,-1.96,9.15,13.23,2.2,3.33,-4.11,3.15,-13.49,-6.2,3.28,-2.79,4.0,-7.25,-10.58,6.6,-3.99,5.19,-2.2,2.43,-5.68,5.45]
    Δφ_Co_60_pure_15000 = [-2.22,-10.11,1.44,-11.69,-2.81,3.81,-5.59,2.96,1.78,-2.82,-5.12,-2.92,-3.87,0.76,4.01,-0.29,9.1,12.52,1.48,2.65,-4.47,2.06,-14.09,-7.62,3.26,-1.73,3.22,-7.62,-9.48,7.2,-1.83,6.47,-2.35,1.32,-6.32,5.93]
    #define constants
    mean = np.mean(Co60_10000_events)
    variance = np.var(Co60_10000_events)
    sigma = np.sqrt(variance)
    mean_of_events = np.mean(10000)
    
    LSL = min(Co60_10000_events)
    USL = max(Co60_10000_events)
    num_bins = 80
    # calculate the z-transform
    #z1 = ( LSL - mean ) / sigma
    #z2 = ( USL - mean ) / sigma

    #Process capability
    Cp = (USL - LSL)/ (6*variance)
    #Print info table
    #info_table = [['delta_phi_upper'],['Entries' , len(real_data)], ['Mean', mean] ,['Variance', variance] ,['Sigma', sigma]]
    #print(tabulate(info_table,headers = 'firstrow', tablefmt = 'fancy_grid'))


    fig, ax = plt.subplots(1,1,figsize=(10,7))
    plt.style.use('seaborn-bright')
    #ax.plot(x_all,y_all, label='Gaussian Distribution', linewidth = 2, color = 'r')
    #ax.legend(title ='delta_phi_upper')
    n, bins, patches = ax.hist(Co60_10000_events, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
    
    
    # add a 'best fit' line
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mean))**2))
    ax.plot(bins, y, '--',label='best fit', color='red')
    ax.set_xlim([-200,200])
    #ax.set_ylim([0,0.45])
    ax.set_xlabel('Δφ[deg]')
    ax.set_ylabel('Events / 1.5deg')
    ax.set_title('Co60 Source (Energy = 1173 keV, 1332 keV) φ Accuracy from upper SiPM')
    plt.annotate("Entries = {}".format(round(mean_of_events)), xy=(0, 1), xytext=(12, -12), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Mean = {0:.4g}".format(mean), xy=(0, 1), xytext=(12, -24), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Variance = {0:.4g}".format(variance), xy=(0, 1), xytext=(12, -36), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Sigma = {0:.4g}".format(sigma), xy=(0, 1), xytext=(12, -48), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.annotate("Process Capability = {0:.4g}".format(Cp), xy=(0, 1), xytext=(12, -60), va='top',xycoords='axes fraction', textcoords='offset points')
    plt.show()


plot()

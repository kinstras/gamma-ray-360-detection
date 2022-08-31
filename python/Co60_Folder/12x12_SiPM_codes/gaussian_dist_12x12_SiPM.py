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

    Co60_3500_events_btin_1150_1350kev = [-17.77,-5,-2.59,28.23,-47.84,4.76,25.37,13.79,17.29,26.75,-30.32, -158.14,5.15,-46.44,13.32,46.69,-30.31,31.08,24.22,-13.38,15.98,45.31,27.22,17.65,-17.77,-10.83,20.16,27.96,17.81,31.04,6.51,-10.28,-5.94,-2.13,2.68,-2.07]
    Co60_7000_events_btin_550_1350kev = [0.35,1.71,1.9,17.85,-12.38,12.36,12.04,7.96,11.56,22.53,-3.99,17.62,4.6,-18.46,3.79,10.41,-10.24,21.83,4.86,-19.37,16.37,30,5.06,12.04,-6.49,-10,15.13,9.14,15.21,0.28,-7.43,-3.52,-10.21,4.89,-0.29]
    Co60_10000_events_btin_0_1350kev = [1.89, 0.4, -1.42, 16.95, -16.22, 11.4, 9.28, 8.78, 6.03, 18.04, -3.42, 19.87, 3.42, -17.34, 4.94, 8.9, -10.92, 20.32, 4.79, -23.29, 15.54, 30.17, 1.93, 10.97, -5.62, -2.53, 14.55, -7.28, 9.23, 15.6, 2.82, -4.43, -2.34, -7.06, 2.34, 2.82]
    #define constants
    mean = np.mean(Co60_10000_events_btin_0_1350kev)
    variance = np.var(Co60_10000_events_btin_0_1350kev)
    sigma = np.sqrt(variance)
    mean_of_events = np.mean(7000)
    
    LSL = min(Co60_10000_events_btin_0_1350kev)
    USL = max(Co60_10000_events_btin_0_1350kev)
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
    n, bins, patches = ax.hist(Co60_10000_events_btin_0_1350kev, bins = num_bins, range=(-200,200),density=True,  edgecolor="blue",color='white')
    
    
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

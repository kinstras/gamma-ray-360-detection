#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import array
import seaborn as sb
from pathlib import Path
import sys
import os
import shutil
import plotly.express as px
import config
from co57_isotope import Co57_isotope
from cs137_isotope import Cs137_isotope
from co60_isotope import Co60_isotope




def final_resolution_plot_for_each_arrangement():
    fig,  ax = plt.subplots()
    #define x and y values
    x_number_of_events = [3500, 5000, 10000, 15000]
    y_resolution_Co57_8x8_SiPM = [1.632, 0.857, 0.89, 0.693]
    y_resolution_Cs137_8x8_SiPM = [11.84, 9.71, 6.4, 5.054]
    y_resolution_Co60_8x8_SiPM = [17.91,15.44, 11.24, 8.847]

    y_resolution_Co57_12x12_SiPM = [1.364, 1.426, 1.04, 0.931]
    y_resolution_Cs137_12x12_SiPM = [10.12, 10.06, 6.257, 5.056]
    y_resolution_Co60_12x12_SiPM = [23.53, 18.77, 11.51, 9.57]
    # Define each line separately
    line1, = plt.plot(x_number_of_events, y_resolution_Co57_12x12_SiPM, c='r',marker="v",ls='-', label='$Co^{57}$')
    line2, = plt.plot(x_number_of_events, y_resolution_Cs137_12x12_SiPM, c='g',marker=(8,2,0),ls='--', label='$Cs^{137}$')
    line3, = plt.plot(x_number_of_events, y_resolution_Co60_12x12_SiPM, c='m',marker="o",ls='--',fillstyle='none', label='$Co^{60}$')

    ax.set_yticks(np.arange(0,max(y_resolution_Co60_8x8_SiPM)+3,2.5))

    plt.legend(loc='upper right')
    #plt.title("Γράφημα resolution της διάταξης \n 8x8 ανιχνευτών πιριτίου ")
    #plt.title("Γράφημα resolution της διάταξης \n 12x12 ανιχνευτών πιριτίου ")
    plt.xlabel("Αριθμός γ-φωτονίων Φωτοκορυφής")
    plt.ylabel("Ακρίβεια αζιμουθιακής γωνίας [$^{\circ}$]")
    #plt.title("Γράφημα resolution της διάταξης 12x12 ανιχνευτών πιριτίου ")
    plt.show()

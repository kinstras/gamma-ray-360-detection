import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tensorflow as tf
from sklearn import model_selection
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dense
from sklearn.preprocessing import StandardScaler
############## proergasia
angle = []
for i in range(0,100,10):
    for j in range(0,1000):
        angle.append(i)
        
df_angle = pd.DataFrame(angle)   


#file_name= input("Enter file name: \t")
file_name = "merge.txt"
used_columns = list(range(144))
input_data = np.genfromtxt(file_name , delimiter = ' ' , filling_values = 0, usecols=(used_columns))
df_raw_data = pd.DataFrame(input_data)  
#print(df_raw_data)
final_df = pd.concat([df_raw_data, df_angle], axis=1, join='inner')


x = final_df.iloc[:,0:-1]
y = final_df.iloc[:,-1]

x_train, x_test, y_train, y_test = model_selection.train_test_split(x, y, test_size=0.1, random_state=4)
########## normalize data
xnorm = StandardScaler();
ynorm = StandardScaler();
x_train=xnorm.fit_transform(x_train)
x_test=xnorm.transform(x_test)
y_train=ynorm.fit_transform(np.array(y_train).reshape(-1,1))
y_test=ynorm.transform(np.array(y_test).reshape(-1,1))
















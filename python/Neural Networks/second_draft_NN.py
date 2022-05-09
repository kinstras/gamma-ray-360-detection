# -*- coding: utf-8 -*-
"""
Created on Mon May  9 23:05:30 2022

@author: kostm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tensorflow as tf
from sklearn import model_selection
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dense
from sklearn.preprocessing import StandardScaler
import math
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
y = y / 10   #to 1 -> 10 moires kok

x_train, x_test, y_train, y_test = model_selection.train_test_split(x, y, test_size=0.1, random_state=4)

#convering x_train %% x_test 0-1
x_train = x_train/ 143
x_test = x_test/ 143

print("First Label before conversion:")
print(y_train[2])
# Converting labels to one-hot encoded vectors
y_train = tf.keras.utils.to_categorical(y_train)
y_test = tf.keras.utils.to_categorical(y_test)

print("First Label after conversion:")
print(y_train[2])
####### neural net
model = tf.keras.Sequential([
    tf.keras.layers.Dense(512, input_shape=(144,), activation='relu'),
    tf.keras.layers.Dense(256, activation='relu'),
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dense(64, activation='relu'),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dense(16, activation='relu'),
    tf.keras.layers.Dense(10, activation='softmax')
])

model.compile(
  loss = 'categorical_crossentropy',
  optimizer = 'adam',
  metrics = ['accuracy']
)
#print(model.summary())
history = model.fit(
  x = x_train,
  y = y_train,
  epochs = 30
)

# # Showing plot for loss
# plt.plot(history.history['loss'])
# plt.xlabel('epochs')
# plt.legend(['loss'])
# plt.show()
 
# # Showing plot for accuracy
# plt.plot(history.history['accuracy'], color='orange')
# plt.xlabel('epochs')
# plt.legend(['accuracy'])
# plt.show()

# Call evaluate to find the accuracy on test images
test_loss, test_accuracy = model.evaluate(
  x = x_test, 
  y = y_test
)
 
print("Test Loss: %.4f"%test_loss)
print("Test Accuracy: %.4f"%test_accuracy)




predicted_probabilities = model.predict(x_test)
predicted_classes = tf.argmax(predicted_probabilities, axis=-1).numpy()

def predict():
    index=2
     
    
    # Printing Probabilities
    print("Probabilities predicted for image at index", index)
    print(predicted_probabilities[index])
     
    print()
     
    # Printing Predicted Class
    print("Probabilities class for image at index", index)
    print(predicted_classes[index]*10 , " moires")


predict()

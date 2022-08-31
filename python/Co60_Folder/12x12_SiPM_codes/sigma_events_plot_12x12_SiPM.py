# Importing packages
import matplotlib.pyplot as plt

# Define x and y values
x_realistic_8x8_SiPM = [7000,10000,15000,20000]
y_realistic_8x8_SiPM = [11.75,8.328,5.861,5.949]

x_realistic_12x12_SiPM = [3500,7000,10000]
y_realistic_12x12_SiPM = [35.28,11.7,11.3]
# Plot a simple line chart without any feature
plt.plot(x, y)
plt.title('Relation Sigma - events')
plt.xlabel('Number of events')
plt.ylabel('Sigma')
plt.show()

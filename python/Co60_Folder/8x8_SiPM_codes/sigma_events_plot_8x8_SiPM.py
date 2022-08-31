# Importing packages
import matplotlib.pyplot as plt

# Define x and y values
x = [7000,10000,15000,20000]
y = [11.75,8.328,5.861,5.949]

# Plot a simple line chart without any feature
plt.plot(x, y)
plt.title('Relation Sigma - events')
plt.xlabel('Number of events')
plt.ylabel('Sigma')
plt.show()

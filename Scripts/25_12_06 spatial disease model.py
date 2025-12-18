import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random
def decision(probability):
    return random.random() < probability

hosts = 100 # number of hosts
fieldsize = 10 # side length of field
timesteps = 500

infected = np.zeros(hosts)

rng = np.random.default_rng(12345)
hostCoords = rng.random((hosts, 2)) * fieldsize
hostCoords = pd.DataFrame(hostCoords)

distanceMatrix = np.zeros((hosts, hosts))

for i in range(0,hosts):
    iCoords = hostCoords.loc[i]
    for j in range(0,hosts):
        jCoords = hostCoords.loc[j]

        distanceMatrix[i,j] = np.sqrt((iCoords[0]-jCoords[0])**2 + (iCoords[1]-jCoords[1])**2)

alpha = 0.3

probabilityMatrix = np.exp(-(distanceMatrix / alpha)) # exponential kernel

A = probabilityMatrix.sum() / hosts

probabilityMatrix = probabilityMatrix / A
probabilityMatrix.sum() / hosts # should be 1

# Run model
infected = np.zeros(hosts)
infected[0] = 1
# infected[2] = 1

# infected_test = decision(infected.dot(probabilityMatrix) + infected)
# infected_test = decision(infected_test.dot(probabilityMatrix) + infected_test)
# infected_test

infection_record = np.zeros((timesteps,hosts))
total_infected = np.zeros(timesteps)
for t in range(0, timesteps):
    infection_record[t] = infected
    total_infected[t] = infected.sum()
    infected = decision(infected.dot(probabilityMatrix) + infected)
total_infected

fig, ax0 = plt.subplots(figsize=(10, 7))
ax0.plot(range(0, timesteps), total_infected)
fig.show()


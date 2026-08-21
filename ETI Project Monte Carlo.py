import numpy as np
import matplotlib.pyplot as plt
import time
import math

# start calculation of simulation time

startTime = time.perf_counter()

# physical parameters

# time data

ThalfLife = 12.33 * 365 # years -> days
rodTime = 365 # days

# tritium production data

TProduced = 400 # grams
TRodSimSeg = 0.001 
TmolarMass = 3.016049 # grams/mol

# permeation data

perm = 0.000436 # grams/day
permVariance = 9.17 * 10**(-5) # +- grams/day

# detector data

detectorEfficiency = 0.466
detectorEfficiencyStDev = 0.002 / 2

# simulation parameters

testSegments = math.floor(TProduced / (TRodSimSeg * TmolarMass)) # number of test segments, translates roughly to 0.001 mol per segment
permSeg = perm * testSegments / TProduced # converts permeation rate to segments
permVarSeg = permVariance * testSegments / TProduced # converts permeation variance to segments
sims = 1000

# monte carlo implementation, thank you ETI staff!!!

# determine decay probability

def decayProb(halfLife, dt):
    # calculates decay probability over timestep

    decay_rate = np.log(2) / halfLife
    return 1 - np.exp(-decay_rate * dt)

# monte carlo function

def monteCarlo(N0, halfLife, t, dt):

    # simulates the decay of N0 segments over time t with timesteps dt
   
    time_grid = np.arange(0, t + dt, dt)
    N_vals = np.zeros_like(time_grid, dtype=int)

    N = N0
    N_vals[0] = N
    
    lambdaD = decayProb(halfLife, dt)

    for t in range(1, len(time_grid)):
        if N > 0:
            # generates N random floats between 0 and 1 then checks against decay probability, also calculates the amount permeated through target walls

            pDecayed = np.random.rand(N) < lambdaD
            perm = permSeg + np.random.uniform(-permVarSeg, permVarSeg)

            # subtracts decayed segments and permeated segments from total remaning tritium

            N = math.floor(N - np.sum(pDecayed) - perm)

        N_vals[t] = N

    return N_vals, time_grid

plt.figure()

# run loss sims, running loss on mol rather than atoms, cause thats way too many atoms

finalValues = [0] * sims

for i in range(sims):
    vals, t = monteCarlo(testSegments, ThalfLife, rodTime, 1) 

    plt.plot(t, TProduced - (vals * TProduced / testSegments))    
    
    finalValues[i] = (testSegments - vals[rodTime]) * TProduced / testSegments # convert lost tritium from segments to grams
 
# generate statistical data

avg = sum(finalValues) / sims# find average
stdev = np.std(finalValues)

endTime = time.perf_counter()

# report results and simulation time

print("\nAverage Lost Tritium: " + str(avg) + " g\nStandard Deviation: " + str(stdev) + " g\nUncertainty (95% Confidence Interval): +- " + str(stdev * 2) + " g\nSimulation Time: " + str(endTime - startTime) + " s\n")

# draw plot

plt.title("Tritium Loss Over a Year (Monte Carlo Simulation)")
plt.xlabel("Time (d)")
plt.ylabel("Lost Tritium (g)")
plt.grid(True)

plt.show()





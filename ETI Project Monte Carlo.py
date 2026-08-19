import numpy as np
import matplotlib.pyplot as plt
import time
import math

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
sims = 100

# monte carlo implementation, thank you ETI staff!!!

# calculation of decay probabilities

def decayProb(halfLife, dt):
    """Calculates decay probability p over time step dt."""
    decayRate = np.log(2) / halfLife
    return 1 - np.exp(-decayRate * dt)

# monte carlo function

def decayProb(halfLife, dt):
    """Calculates decay probability p over time step dt."""
    decay_rate = np.log(2) / halfLife
    return 1 - np.exp(-decay_rate * dt)

def monteCarlo(N_0, halfLife, t_tot, dt):
    """Simulates decay of N_0 particles over total time t_tot in steps of dt."""
   
    time_grid = np.arange(0, t_tot + dt, dt)
    N_vals = np.zeros_like(time_grid, dtype=int)

    N = N_0
    N_vals[0] = N
    
    lambdaD = decayProb(halfLife, dt)

    for t in range(1, len(time_grid)):
        if N > 0:
            # Generate N random floats between 0 and 1; check against decay probability
            pDecayed = np.random.rand(N) < lambdaD
            perm = permSeg + np.random.uniform(-permVarSeg, permVarSeg)

            N = math.floor(N - np.sum(pDecayed) - perm)

        N_vals[t] = N

    return N_vals, time_grid

noDetectAvg = 0
#detectAvg = 0

plt.figure()

# run loss sims, running loss on mol rather than atoms, cause thats way too many atoms

for i in range(sims):
    vals, t = monteCarlo(testSegments, ThalfLife, rodTime, 1) 

    plt.plot(t, TProduced - (vals * TProduced / testSegments))    
    
    noDetectAvg = noDetectAvg + vals[rodTime]
 
    #detectAvg = detectAvg + detectVals[rodTime]

# generate statistical data

noDetectAvg = noDetectAvg * TProduced / testSegments # convert from segments to grams
noDetectAvg = TProduced - (noDetectAvg / sims) # find amount lost
noDetectStdev = np.sqrt(noDetectAvg)

#detectAvg = detectAvg / sims
#detectStdev = np.sqrt(detectAvg)

# report results

print("\nAverage Lost Tritium (No Detector): " + str(noDetectAvg) + "g\nStandard Deviation (No Detector): " + str(noDetectStdev) + "g\nUncertainty (No Detector, 95% Confidence Interval): +-" + str(noDetectStdev * 2) + "g\n")
# print("\n\nAverage Lost Tritium (With Detector): " + str(detectAvg) + "\nStandard Deviation (With Detector): " + str(detectStdev))

# determine hiding space between detected and simulated, propagating uncertainty according to sqrt(sigmaA^2 + sigmaB^2)

# print("\n\nDifference Between Detected and Simulated: " + str(noDetectAvg - detectAvg) + "\nStandard Deviation (With Detector): " + str(np.sqrt(detectStdev**2 + noDetectStdev**2)))
# print("\n\nEfficiency Scaled Estimate: " + str(detectAvg / detectorEfficiency) + "\nStandard Deviation (Scaled): " + str(np.sqrt(detectStdev**2 + detectorEfficiencyStDev**2)))
# print("\n\nDifference Between Scaled Detected and Simulated: " + str(noDetectAvg - (detectAvg / detectorEfficiency)) + "\nStandard Deviation (Scaled): " + str(np.sqrt(detectStdev**2 + noDetectStdev**2 + noDetectStdev**2)))
# print("\nAverage Unaccounted For Tritium: " + str(noDetectAvg - (detectAvg / detectorEfficiency)) + " +- " + str(np.sqrt(detectStdev**2 + noDetectStdev**2 + noDetectStdev**2) * 2) + " (95 percent confidence)\n")

# draw plots

plt.title("Loss Simulations (No Detector)")
plt.xlabel("Time (d)")
plt.ylabel("Lost Tritium (g)")

plt.show()





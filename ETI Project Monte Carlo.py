import numpy as np
import matplotlib.pyplot as plt
import time

# material parameter design

ThalfLife = 12.33 # years
hWTemp = 
fuelTemp = 
airTemp = 
reacPress =
airPress = 
rodRadius = 
airLen = 
hWLen = 
fuelLen = 
TRodSimMol = 0.001 # mol
TmolarMass = 3.016049 # g/mol
rodTime = 1 # year
porosity =
permeability =
zirThickness = 
detectorEfficiency = 46.6
detectorEfficiencyStDev = 0.2 / 2

sims = 30

# monte carlo implementation, thank you Gracie!!!

def decayProb(halfLife, dt):
    """Calculates decay probability p over time step dt."""
    decayRate = np.log(2) / halfLife
    return 1 - np.exp(-decayRate * dt)

def permProb():


def monteCarlo(N0, halfLife, tTot, dt):
    """Simulates decay of N_0 particles over total time t_tot in steps of dt."""
    timeGrid = np.arange(0, tTot + dt, dt)
    NVals = np.zeros_like(timeGrid, dtype=int)
    NValsDetect = NVals

    N = N0
    NVals[0] = N
    NValsDetect[0] = N
    lambdaD = decayProb(halfLife, dt)
    lambdaP = permProb()
    lambdaDet = detectorEfficiency

    for t in range(1, len(timeGrid)):
        if N > 0:
            # Generate N random floats between 0 and 1; check against decay probability
            pDecayed = np.random.rand(N) < lambdaD
            pPerm = np.random.rand(N) < lambdaP
            pDetect = np.random.rand(N) < lambdaDet

            N = N - np.sum(pDecayed) - np.sum(pPerm)
        NVals[t] = N
        NValsDetect[t] = N - np.sum(pDetect)

    return NVals, timeGrid

noDetectAvg = 0
detectAvg = 0

fig, (noDetect, detect) = plt.subplots(1, 2, figsize=(10, 4))

# run  loss sims

for i in range(sims):
    vals, detectVals, t = monteCarlo()
    vals = vals * molarMass / (6.02 * 10**23) # convert to grams  
    vals = vals * multFactor # convert to total grams in all rods
    noDetect.plot(t, vals)
    noDetectAvg = noDetectAvg + vals[TRodSimAtoms]

    detectVals = detectVals * molarMass / (6.02 * 10**23) # convert to grams  
    detectVals = detectVals * multFactor # convert to total grams in all rods
    detectAvg = detectAvg + detectVals[TRodSimAtoms]

# generate statistical data

noDetectAvg = noDetectAvg / sims
noDetectStdev = np.sqrt(noDetectAvg)

detectAvg = detectAvg / sims
detectStdev = np.sqrt(detectAvg)

# report results

print("Average Lost Tritium (No Detector): " + str(noDetectAvg) + "\nStandard Deviation (No Detector): " + str(noDetectStdev))
print("\n\nAverage Lost Tritium (With Detector): " + str(detectAvg) + "\nStandard Deviation (With Detector): " + str(detectStdev))

# determine hiding space between detected and simulated, propagating uncertainty according to sqrt(sigmaA^2 + sigmaB^2)

print("\n\nDifference Between Detected and Simulated: " + str(noDetectAvg - detectAvg) + "\nStandard Deviation (With Detector): " + str(np.sqrt(detectStdev^2 + noDetectStdev^2)))
print("\n\nEfficiency Scaled Estimate: " + str(detectAvg / detectorEfficiency) + "\nStandard Deviation (Scaled): " + str(np.sqrt(detectStdev^2 + detectorEfficiencyStDev^2)))
print("\n\nDifference Between Scaled Detected and Simulated: " + str(noDetectAvg - (detectAvg / detectorEfficiency)) + "\nStandard Deviation (Scaled): " + str(np.sqrt(detectStdev^2 + noDetectStdev^2 + noDetectStdev^2)))
print("\nAverage Unaccounted For Tritium: " + str(noDetectAvg - (detectAvg / detectorEfficiency)) + " +- " + str(np.sqrt(detectStdev^2 + noDetectStdev^2 + noDetectStdev^2) * 2) + " (95 percent confidence)")

# draw plots

noDetect.title("Loss Simulations (No Detector)")
noDetect.xlabel("Time")
noDetect.ylabel("Remaining Tritium")

noDetect.title("Loss Simulations (Detector)")
noDetect.xlabel("Time")
noDetect.ylabel("Remaining Tritium")

plt.tight_layout
plt.show





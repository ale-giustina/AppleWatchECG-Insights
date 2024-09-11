import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 

ecgfiles = os.listdir('ECG')

def peak_detector(array, threshold_min=None, threshold_max=None, inverted=False):
    peaks = []
    
    if inverted:
        
        for i in range(1, len(array)-1):
            if array[i] <= array[i-1] and array[i] <= array[i+1]:
                if threshold_max is not None and threshold_min is not None:
                    if array[i] > threshold_min and array[i] < threshold_max:
                        peaks.append(i)
                else:
                    peaks.append(i)
    
    else:
        
        for i in range(1, len(array)-1):
            if array[i] >= array[i-1] and array[i] >= array[i+1]:
                if threshold_max is not None and threshold_min is not None:
                    if array[i] > threshold_min and array[i] < threshold_max:
                        peaks.append(i)
                else:
                    peaks.append(i)
    return peaks

def derivative(array):
    deriv = np.zeros(len(array)-1)
    for i in range(len(array)-1):
        deriv[i] = array[i+1] - array[i]
    return deriv

def check_dist(point, peaks, dist=50):
    for peak in peaks:
        if abs(point-peak) < dist:
            return False
    return True


df = pd.read_csv('ECG/' + ecgfiles[0], header=None, skiprows=12)

df = df [0]

array = np.array(df)


avg_array = array


peaks_deriv_r = peak_detector(derivative(avg_array),-1000,-100,inverted=True)

min_dist = 50

for i in range(1, len(peaks_deriv_r)-1):
    if peaks_deriv_r[i] - peaks_deriv_r[i-1] < min_dist:
        peaks_deriv_r[i] = 0
if peaks_deriv_r[-1] - peaks_deriv_r[-2] < min_dist:
    peaks_deriv_r[-1] = 0
peaks_deriv_r = [i for i in peaks_deriv_r if i != 0]

#find r peaks based on the derivative
r_peak=[]

for indx, l in enumerate(peaks_deriv_r):
    
    prev_val = avg_array[l]  
    
    while l > 0:  
        l -= 1  
        if avg_array[l] < prev_val: 
            r_peak.append(l+1)  
            break
        prev_val = avg_array[l]

r_peak=list(dict.fromkeys(r_peak))

# #find T peaks based on derivatives
# T_peak=[]
# peaks_deriv_T = peak_detector(derivative(avg_array),-290,-50, inverted=True)

# # max_dist=0
# T_peak = []

# for indx, l in enumerate(peaks_deriv_T):
#     max_dist = 0
#     prev_val = avg_array[l]  
#     while l > 0:  
#         max_dist += 1
#         l -= 1  
#         if avg_array[l] < prev_val and check_dist(l, r_peak, 20): 
#             T_peak.append(l+1)  
#             break
        
#         prev_val = avg_array[l]  


#         if max_dist > 10:
#             break

#find T peaks based on max value between two R peaks
T_peak = []

for indx, l in enumerate(r_peak[:-1]):
    start_indx = l+15
    end_indx = r_peak[indx+1]-15
    T_peak.append(start_indx + np.argmax(avg_array[start_indx:end_indx]))

fig, ax = plt.subplots(2,1,sharex=True)
ax[0].plot(avg_array)
ax[0].plot(r_peak, [avg_array[i] for i in r_peak], 'go')
ax[0].plot(T_peak, [avg_array[i] for i in T_peak], 'ro')
ax[0].grid()

ax[1].plot(derivative(avg_array))
ax[1].plot(peaks_deriv_r, [derivative(avg_array)[i] for i in peaks_deriv_r], 'ro')
#ax[1].plot(peaks_deriv_T, [derivative(avg_array)[i] for i in peaks_deriv_T], 'bo')

ax[1].grid()


plt.show()
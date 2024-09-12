import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Formatter
import pandas as pd
import os

def find_square(area):
    width = int((np.sqrt(area))+0.5)
    for i in range(width, 0, -1):
        if area % i == 0:
            return i, area//i

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



filenum = 5
filepath="ECG2"

showentire = False
showcuts = False

saveentire = True
savecuts = True

ecgfiles = os.listdir(filepath)
df = pd.read_csv(f'{filepath}/' + ecgfiles[filenum], header=None, skiprows=12)

df = df [0]
array = np.array(df)

der_array = derivative(array)


#find R peaks based on the derivative
peaks_deriv_r = peak_detector(der_array,-1000,-80,inverted=True)


min_dist = 50

for i in range(1, len(peaks_deriv_r)-1):

    if peaks_deriv_r[i] - peaks_deriv_r[i-1] < min_dist:

        peaks_deriv_r[i] = 0

if peaks_deriv_r[-1] - peaks_deriv_r[-2] < min_dist:

    peaks_deriv_r[-1] = 0

peaks_deriv_r = [i for i in peaks_deriv_r if i != 0]


#find true r peaks based on the derivative
R_peaks=[]


for indx, l in enumerate(peaks_deriv_r):

    prev_val = array[l]

    while l > 0:

        l -= 1

        if array[l] < prev_val:

            R_peaks.append(l+1)
            break
        
        prev_val = array[l]

R_peaks=list(dict.fromkeys(R_peaks))



#find T peaks based on the min value between two R peaks in the derivative
T_peak_deriv = []


for indx, l in enumerate(peaks_deriv_r[:-1]):

    start_indx = l+50

    end_indx = peaks_deriv_r[indx+1]-(peaks_deriv_r[indx+1]-start_indx)//2

    T_peak_deriv.append(start_indx + np.argmin(der_array[start_indx:end_indx]))



#find true T peaks position based on the derivative by going back
T_peaks=[]
T_peaks_end = []

for indx, l in enumerate(T_peak_deriv):

    strt_in = l

    prev_val = array[strt_in]

    #find true T peak

    while strt_in > 0:

        strt_in -= 1

        if array[strt_in] < prev_val:

            T_peaks.append(strt_in+1)
            break
        
        prev_val = array[strt_in]
    
    #find end peak with normal graph

    # while l > 0:

    #     l += 1

    #     if array[l] > prev_val:

    #         T_peaks_end.append(l-1)
    #         break
        
    #     prev_val = array[l]

    #find end peak with derivative

    prev_val = der_array[l]

    num_of_t = 0

    while l > 0:

        l += 1

        if der_array[l] < prev_val or der_array[l] > -5:
            if num_of_t > 3:
                T_peaks_end.append(l+6)
                break
            else:
                num_of_t+=1
        
        prev_val = der_array[l]

#make list unique
T_peaks=list(dict.fromkeys(T_peaks))


#Find Q and S peaks based on R peaks
Q_peaks = []
Q_start_peaks = []
S_peaks = []
S_end_peaks = []

for indx, R_peak in enumerate(R_peaks):

    l = R_peak

    Q_peak_prev = array[l]
    S_peak_prev = array[l]

    #find the actual peak

    while l > 0:

        l-=1

        if array[l] > Q_peak_prev:

            Q_peaks.append(l+1)
            break
        
        Q_peak_prev = array[l]

    #find the start of the peak

    dist = 0
    while l > 0:

        l-=1

        if array[l] < Q_peak_prev:

            Q_start_peaks.append(l+1)
            break
        
        #if start of Q peak is blurred with P wave

        if dist>15 and array[l] - Q_peak_prev < 0.5:
            
            Q_start_peaks.append(l+1)
            break

        Q_peak_prev = array[l]
        dist+=1

    l = R_peak

    #find the actual peak

    while l < len(array)-1:

        l+=1

        if array[l] > S_peak_prev:

            S_peaks.append(l-1)
            break
        
        S_peak_prev = array[l]

    #find the end of the peak

    while l < len(array)-1:

        l+=1

        if array[l] < S_peak_prev:

            S_end_peaks.append(l-1)
            break
        
        S_peak_prev = array[l]

#group toghether the R T Q S peaks into each beat
detected_peaks = []

for index, Q_start in enumerate(Q_start_peaks):
    
    try:
        detected_peaks.append([Q_start, Q_peaks[index], R_peaks[index], S_peaks[index], S_end_peaks[index], T_peaks[index], T_peaks_end[index]])
    
    except:
        pass



fig, ax = plt.subplots(2,1,sharex=True, figsize=(100,30))

ax[0].plot(array)
ax[0].plot(R_peaks, [array[i] for i in R_peaks], 'go', label='R peaks')
ax[0].plot(T_peaks, [array[i] for i in T_peaks], 'ro', label='T peaks')
ax[0].plot(T_peaks_end, [array[i] for i in T_peaks_end], 'o', label='T peaks', color="black")
ax[0].plot(Q_peaks, [array[i] for i in Q_peaks], 'co', label='Q peaks')
ax[0].plot(S_peaks, [array[i] for i in S_peaks], 'bo', label='S peaks')
ax[0].plot(S_end_peaks, [array[i] for i in S_end_peaks], 'o', label='S end peaks', color="gray")
ax[0].plot(Q_start_peaks, [array[i] for i in Q_start_peaks], 'o', label='Q start peaks', color="orange")
ax[0].legend()
ax[0].grid()


ax[1].plot(der_array)
ax[1].plot(peaks_deriv_r, [der_array[i] for i in peaks_deriv_r], 'ro', label='R peaks')
ax[1].plot(T_peak_deriv, [der_array[i] for i in T_peak_deriv], 'bo', label='T peaks')
ax[1].legend()
ax[1].grid()

if saveentire:
    plt.savefig('Examples/' + ecgfiles[filenum].split('.')[0] + '.png')
    print('saved: '+ ecgfiles[filenum].split('.')[0] + '.png')

if showentire:
    plt.show()

plt.close()


if len(detected_peaks) <= 49:
    wid, hei = (7,7)
if len(detected_peaks) <= 36:
    wid, hei = (6,6)
if len(detected_peaks) <= 25:
    wid, hei = (5,5)
if len(detected_peaks) <= 16:
    wid, hei = (4,4)

inx, iny = 0, 0

fig, ax = plt.subplots(wid, hei, figsize=(50,30))

for index, peak in enumerate(detected_peaks):
    
    ax[inx, iny].plot(array[peak[0]:peak[6]])
    ax[inx, iny].plot(0, array[peak[0]], 'o', label='Start Q peak', color="orange")
    ax[inx, iny].plot(peak[1]-peak[0], array[peak[1]], 'co', label='Q peak')
    ax[inx, iny].plot(peak[2]-peak[0], array[peak[2]], 'go', label='R peak')
    ax[inx, iny].plot(peak[3]-peak[0], array[peak[3]], 'bo', label='S peak')
    ax[inx, iny].plot(peak[4]-peak[0], array[peak[4]], 'o', label='S end peak', color="gray")
    ax[inx, iny].plot(peak[5]-peak[0], array[peak[5]], 'ro', label='T peak')
    ax[inx, iny].plot(peak[6]-peak[0], array[peak[6]], 'o', label='T peak end', color="black")
    ax[inx, iny].grid()

    iny+=1

    if iny == hei:
        iny = 0
        inx+=1

fig.legend(['microvolt','Q start peak','Q peak','R peak','S peak','S end peak','T peak','T peak end'], loc='upper right')

if showcuts:
    plt.show()
if savecuts:
    plt.savefig('Examples2/' + ecgfiles[filenum].split('.')[0] + '.png')
    print('saved: '+ ecgfiles[filenum].split('.')[0] + '.png')
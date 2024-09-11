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

filenum = 5

df = pd.read_csv('ECG/' + ecgfiles[filenum], header=None, skiprows=12)


df = df [0]


array = np.array(df)

der_array = derivative(array)


peaks_deriv_r = peak_detector(der_array,-1000,-80,inverted=True)


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

    prev_val = array[l]


    while l > 0:

        l -= 1
        if array[l] < prev_val:

            r_peak.append(l+1)

            break
        prev_val = array[l]


r_peak=list(dict.fromkeys(r_peak))


#find T peaks based on the min value between two R peaks in the derivative

T_peak_deriv = []


for indx, l in enumerate(peaks_deriv_r[:-1]):

    start_indx = l+60

    end_indx = peaks_deriv_r[indx+1]-100

    T_peak_deriv.append(start_indx + np.argmin(der_array[start_indx:end_indx]))

#find true T peaks position based on the derivative
T_peak=[]


for indx, l in enumerate(T_peak_deriv):

    prev_val = array[l]

    while l > 0:

        l -= 1
        if array[l] < prev_val:

            T_peak.append(l+1)

            break
        prev_val = array[l]

#make list unique
T_peak=list(dict.fromkeys(T_peak))

fig, ax = plt.subplots(2,1,sharex=True, figsize=(20,10))

ax[0].plot(array)

ax[0].plot(r_peak, [array[i] for i in r_peak], 'go', label='R peaks')

ax[0].plot(T_peak, [array[i] for i in T_peak], 'ro', label='T peaks')

ax[0].legend()
ax[0].grid()


ax[1].plot(der_array)

ax[1].plot(peaks_deriv_r, [der_array[i] for i in peaks_deriv_r], 'ro', label='R peaks')

ax[1].plot(T_peak_deriv, [der_array[i] for i in T_peak_deriv], 'bo', label='T peaks')

ax[1].legend()
ax[1].grid()

plt.savefig('Examples/' + ecgfiles[filenum].split('.')[0] + '.png')

plt.show()
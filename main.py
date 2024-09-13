import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Formatter
import pandas as pd
import os
import seaborn as sns
from matplotlib.patches import Rectangle
import shutil

plt.switch_backend('agg')

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

def data_extractor(peaks, sample_ratem, remove_outliers=0):

    #peaks structure: [[Q_start, Q, R, S, S_end, T, T_end],[...],...]
    #                   0        1  2  3  4      5  6  
    Sinus_RR = []
    Frequency = 0
    interval_ST = []
    interval_QT = []
    interval_QRS = []

    for i in range(len(peaks)-1):
        
        if i != len(peaks)-1:
            Sinus_RR.append((peaks[i+1][2]-peaks[i][2])/sample_rate)
        else:
            Sinus_RR.append((peaks[i][2]-peaks[i-1][2])/sample_rate)

        interval_ST.append((peaks[i][6]-peaks[i][4])/sample_rate)

        interval_QT.append((peaks[i][6]-peaks[i][0])/sample_rate)

        interval_QRS.append((peaks[i][4]-peaks[i][0])/sample_rate)
    
    Frequency = 60/(np.mean(Sinus_RR))

    if remove_outliers != 0:
        sinus_rr_mean = np.mean(Sinus_RR)
        sinus_rr_std = np.std(Sinus_RR)
        interval_st_mean = np.mean(interval_ST)
        interval_st_std = np.std(interval_ST)
        interval_qt_mean = np.mean(interval_QT)
        interval_qt_std = np.std(interval_QT)
        interval_qrs_mean = np.mean(interval_QRS)
        interval_qrs_std = np.std(interval_QRS)
        
        Sinus_RR = [i for i in Sinus_RR if i > sinus_rr_mean-remove_outliers*sinus_rr_std and i < sinus_rr_mean+remove_outliers*sinus_rr_std]
        interval_ST = [i for i in interval_ST if i > interval_st_mean-remove_outliers*interval_st_std and i < interval_st_mean+remove_outliers*interval_st_std]
        interval_QT = [i for i in interval_QT if i > interval_qt_mean-remove_outliers*interval_qt_std and i < interval_qt_mean+remove_outliers*interval_qt_std]
        interval_QRS = [i for i in interval_QRS if i > interval_qrs_mean-remove_outliers*interval_qrs_std and i < interval_qrs_mean+remove_outliers*interval_qrs_std]
    
    return Frequency, Sinus_RR, interval_ST, interval_QT, interval_QRS



#SETTINGS

#filename or index of the file in the folder
filenum = "ecg_2024-09-12_4.csv"

#folder path
filepath="ECG2"

#show plots interactively
showentire = False
showcuts = False
showanalysis = False
showpoincare = False

#save to folder
save_folder = 'example'
save_svg = False

sample_rate = 512


def analysis(filenum, filepath, showentire=False, showcuts=False, showanalysis=False, showpoincare=False, save_folder=False, save_svg=False, sample_rate=512):

    ecgfiles = os.listdir(filepath)

    if type(filenum) is str:
        try:
            filenum = ecgfiles.index(filenum)
        except:
            raise Exception('File not found: '+filepath+"/"+filenum)

    ECG_filename = ecgfiles[filenum]

    df = pd.read_csv(f'{filepath}/' + ECG_filename, header=None, skiprows=12)

    df = df [0]
    array = np.array(df)

    der_array = derivative(array)

    if save_folder == True:

        if os.path.isdir(ECG_filename.split('.')[0]) == False:
            os.mkdir(ECG_filename.split('.')[0])
        else:
            shutil.rmtree(ECG_filename.split('.')[0])
            os.mkdir(ECG_filename.split('.')[0])

        save_folder = ECG_filename.split('.')[0]

    #find R peaks based on the derivative
    peaks_deriv_r = peak_detector(der_array,-1000,-80,inverted=True)


    min_dist = 50

    for i in range(1, len(peaks_deriv_r)-1):

        if peaks_deriv_r[i] - peaks_deriv_r[i-1] < min_dist:

            peaks_deriv_r[i-1] = 0

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
        
        #find end peak with derivative

        prev_val = der_array[l]

        num_of_t = 0

        while l > 0:

            l += 1

            if der_array[l] < prev_val or der_array[l] > -3:
                if num_of_t > 3:
                    T_peaks_end.append(l+14)
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

    ax[0].set_xticks(np.arange(0,len(array),512/2),)

    #25 ticks per second
    ax[0].set_xticks(np.arange(0,len(array),512/25*5),)
    #convert to seconds
    ax[0].set_xticklabels(np.arange(0,(len(array))/512*1000,1/25*5*1000),)
    ax[1].tick_params(axis='x', rotation=45)
    ax[0].set_xlabel('ms')
    ax[0].set_ylabel('microvolt')
    ax[0].legend()
    ax[0].grid()

    ax[1].plot(der_array)
    ax[1].plot(peaks_deriv_r, [der_array[i] for i in peaks_deriv_r], 'ro', label='R peaks')
    ax[1].plot(T_peak_deriv, [der_array[i] for i in T_peak_deriv], 'bo', label='T peaks')
    ax[1].set_xlabel('ms')
    ax[1].legend()
    ax[1].grid()

    if save_folder:
        plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_entire.png')
        if save_svg:
            plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_entire.svg')
        print('saved: '+ ecgfiles[filenum].split('.')[0] + '_entire.png')

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

    fig.suptitle('Single beats, sample rate: 512Hz')

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

        #25 ticks per second
        ax[inx, iny].set_xticks(np.arange(0,peak[6]-peak[0],512/25),)
        #convert to seconds
        ax[inx, iny].set_xticklabels(np.arange(0,(peak[6]-peak[0])/512*1000,1/25*1000),)
        ax[inx, iny].set_xlabel('ms')
        ax[inx, iny].set_ylabel('microvolt')

        iny+=1

        if iny == hei:
            iny = 0
            inx+=1

    fig.legend(['microvolt','Q start peak','Q peak','R peak','S peak','S end peak','T peak','T peak end'], loc='upper right')

    if showcuts:
        plt.show()
    if save_folder:
        plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_single_peaks.png')
        if save_svg:
            plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_single_peaks.svg')
        print('saved: '+ ecgfiles[filenum].split('.')[0] + '.png')

    plt.close()

    #extract data from peaks and remove outliers above 2 std dev
    Frequency, Sinus_RR, interval_ST, interval_QT, interval_QRS = data_extractor(detected_peaks, sample_rate, remove_outliers=2)

    rr_series = pd.Series(Sinus_RR, name='RR interval')
    st_series = pd.Series(interval_ST, name='ST interval')
    qt_series = pd.Series(interval_QT, name='QT interval')
    qrs_series = pd.Series(interval_QRS, name='QRS interval')

    df = pd.concat([rr_series, st_series, qt_series, qrs_series], axis=1)


    rr_normalrange = [0.6,1.2]

    T_wave = 0.180

    st_normalrange = [(0.08+T_wave),(0.12+T_wave)] #https://www.sciencedirect.com/topics/medicine-and-dentistry/st-segment
    qt_normalrange = [0.350,0.470] #https://www.sciencedirect.com/topics/medicine-and-dentistry/qt-interval
    qrs_normalrange = [0.06,0.12] #https://www.sciencedirect.com/topics/medicine-and-dentistry/qrs-complex

    plt.gca().add_patch(Rectangle((-0.4,rr_normalrange[0]),0.8,rr_normalrange[1]-rr_normalrange[0],linewidth=1,edgecolor='g',facecolor='green', alpha=0.2))
    plt.gca().add_patch(Rectangle((0.6,st_normalrange[0]),0.8,st_normalrange[1]-st_normalrange[0],linewidth=1,edgecolor='g',facecolor='green', alpha=0.2))
    plt.gca().add_patch(Rectangle((1.6,qt_normalrange[0]),0.8,qt_normalrange[1]-qt_normalrange[0],linewidth=1,edgecolor='g',facecolor='green', alpha=0.2))
    plt.gca().add_patch(Rectangle((2.6,qrs_normalrange[0]),0.8,qrs_normalrange[1]-qrs_normalrange[0],linewidth=1,edgecolor='g',facecolor='green', alpha=0.2))
    plt.ylabel('s')

    sns.stripplot(data=df, jitter=0.3, size=3)

    if showanalysis:
        plt.show()
    if save_folder:

        plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_analysis.png')
        if save_svg:
            plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_analysis.svg')
        print('saved: '+ ecgfiles[filenum].split('.')[0] + '_analysis.png')

    plt.close()

    RMSSD = np.sqrt(np.mean(np.square(np.diff(Sinus_RR))))

    #Poincare and variability plot
    fig, ax = plt.subplots(1,1, figsize=(7,7))

    ax.scatter(Sinus_RR[:-1], Sinus_RR[1:], s=15)
    ax.set_xlabel('RR(n)')
    ax.set_ylabel('RR(n+1)')
    ax.set_title('Poincare plot, RMSSD: '+str(np.round(RMSSD*1000,2))+'ms')

    ax.set_ylim(0,1.6)
    ax.set_xlim(0,1.6)

    if showpoincare:
        plt.show()
    if save_folder:
        plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_poincare.png')
        if save_svg:
            plt.savefig(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_poincare.svg')
        print('saved: '+ ecgfiles[filenum].split('.')[0] + '_poincare.png')

    plt.close()

    #create log file

    df = pd.read_csv(f'{filepath}/' + ecgfiles[filenum], header=None)

    #%of datapoints out of norm
    rr_outofnorm = len([i for i in Sinus_RR if i < rr_normalrange[0] or i > rr_normalrange[1]])/len(Sinus_RR)*100
    st_outofnorm = len([i for i in interval_ST if i < st_normalrange[0] or i > st_normalrange[1]])/len(interval_ST)*100
    qt_outofnorm = len([i for i in interval_QT if i < qt_normalrange[0] or i > qt_normalrange[1]])/len(interval_QT)*100
    qrs_outofnorm = len([i for i in interval_QRS if i < qrs_normalrange[0] or i > qrs_normalrange[1]])/len(interval_QRS)*100

    if save_folder:

        with open(f'{save_folder}/' + ecgfiles[filenum].split('.')[0] + '_log.txt', 'w') as f:
            f.write(f'File: {ecgfiles[filenum]}\n')
            f.write(f'{df[0][0]}: {df[1][0]}\n')
            f.write(f'{df[0][1]}: {df[1][1]}\n')
            f.write(f'{df[0][2]}: {df[1][2]}\n')
            f.write(f'{df[0][3]}: {df[1][3]}\n')
            f.write(f'{df[0][4]}: {df[1][4]}\n')
            f.write(f'{df[0][5]}: {df[1][5]}\n')
            f.write(f'{df[0][6]}: {df[1][6]}\n')
            f.write(f'{df[0][7]}: {df[1][7]}\n')
            f.write(f'{df[0][8]}: {df[1][8]}\n')
            f.write(f'{df[0][9]}: {df[1][9]}\n\n')

            f.write(f'avg Frequency: {Frequency}\n\n')
            f.write(f'normal RMSSD interval: 19-107ms\n')
            f.write(f'RMSSD: {np.round(RMSSD*1000,2)}ms\n\n')
            f.write(f'normal RR range: {rr_normalrange[0]*1000}-{rr_normalrange[1]*1000}\n')
            f.write(f'avg Sinus RR: {np.round(np.mean(Sinus_RR)*1000,2)}ms,\nstd dev: {np.round(np.std(Sinus_RR)*1000,2)}ms,\n% of out of norm: {np.round(rr_outofnorm,2)}%\n\n')
            f.write(f'normal ST range: {st_normalrange[0]*1000}-{st_normalrange[1]*1000}\n')
            f.write(f'avg interval ST: {np.round(np.mean(interval_ST)*1000,2)}ms,\nstd dev: {np.round(np.std(interval_ST)*1000,2)}ms,\n% of out of norm: {np.round(st_outofnorm,2)}%\n\n')
            f.write(f'normal QT range: {qt_normalrange[0]*1000}-{qt_normalrange[1]*1000}\n')
            f.write(f'avg interval QT: {np.round(np.mean(interval_QT)*1000,2)}ms,\nstd dev: {np.round(np.std(interval_QT)*1000,2)}ms,\n% of out of norm: {np.round(qt_outofnorm,2)}%\n\n')
            f.write(f'normal QRS range: {qrs_normalrange[0]*1000}-{qrs_normalrange[1]*1000}\n')
            f.write(f'avg interval QRS: {np.round(np.mean(interval_QRS)*1000,2)}ms,\nstd dev: {np.round(np.std(interval_QRS)*1000,2)}ms,\n% of out of norm: {np.round(qrs_outofnorm,2)}%\n\n')

analysis(filenum, filepath, showentire, showcuts, showanalysis, showpoincare, save_folder, save_svg)

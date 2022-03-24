"""
Class definition that implements Heart Rate calculation method.
"""
import numpy as np
import scipy.signal as sig
from math import ceil
import matplotlib.pyplot as plt

class hrMethod():

        # Foresee possibility to include object with all parameters to include
        def __init__(self, fs):
            self.refPeriod = 0.250
            self.THRES = 0.6
            self.fs = fs
            # Parameter to exclude some parts from the analysis
            # Indicate by using indexes
            self.fidVec = np.array([])
            # Force sign of peaks (positive value/negative value)
            self.signForce = np.array([])

            # Constants
            self.MED_SMOOTH_NB_COEFF = round(fs/100)
            self.INT_NB_COEFF = round(7*fs/256)
            self.SEARCH_BACK = 1
            self.MAX_FORCE = np.array([])
            self.MIN_AMP = 0.1
        
        def inputECG(self, ecg):
            """
            A function to rsave ecg data, number of samples and correspondent time.
            
            Parameters
            -----------
            ecg: np.array or list
            
            """
            # Consider numpy array as only acceptable case
            if(type(ecg).__module__ == np.__name__):
                a, b = np.shape(ecg)
                if(a>b):
                    self.num_samples = a
                else:
                    self.num_samples = b
                self.ecg = ecg
            elif(isinstance(ecg, list)):
                self.num_samples = len(ecg)
                self.ecg = np.array(ecg)

            self.time = np.arange(0, self.num_samples)/self.fs

        def filter_band_pass(self):
            """
            Filter the ecg with a sombrero low pass filter.
            """
            b1 = [-7.757327341237223e-5,  -2.357742589814283e-4, -6.689305101192819e-4, -0.001770119249103, \
                    -0.004364327211358, -0.010013251577232, -0.021344241245400, -0.042182820580118, -0.077080889653194, \
                    -0.129740392318591, -0.200064921294891, -0.280328573340852, -0.352139052257134, -0.386867664739069, \
                        -0.351974030208595, -0.223363323458050, 0, 0.286427448595213, 0.574058766243311, \
                    0.788100265785590, 0.867325070584078, 0.788100265785590, 0.574058766243311, 0.286427448595213, 0, \
                    -0.223363323458050, -0.351974030208595, -0.386867664739069, -0.352139052257134, \
                    -0.280328573340852, -0.200064921294891, -0.129740392318591, -0.077080889653194, -0.042182820580118, \
                    -0.021344241245400, -0.010013251577232, -0.004364327211358, -0.001770119249103, -6.689305101192819e-04, \
                    -2.357742589814283e-04, -7.757327341237223e-05]
            # Is resampling of the filter necessary?
            # b1 = sig.resample(b1,len(b1)/(self.fs/250))
            self.filt_ecg = sig.filtfilt(b1,1,self.ecg[:,0])

        def calculateHR(self):
            """
            Main function to calculate HR.
            """

            self.filter_band_pass()

            if (len(np.where(abs(self.filt_ecg)>self.MIN_AMP)[0] )< 0.2):
                # Flat line case
                self.qrs_pos = np.array([])
                self.R_t = np.array([])
                self.R_amp = np.array([])
                self.hrv = np.array([])
                sign = np.array([])
                en_thres = np.array([])
            else:
                diffECG = np.diff(self.filt_ecg)
                sqrECG = np.multiply(diffECG,diffECG)
                intECG = sig.lfilter(np.ones(self.INT_NB_COEFF),1,sqrECG)
                mdfint = sig.medfilt(intECG,kernel_size=self.MED_SMOOTH_NB_COEFF)
                delay = ceil(self.INT_NB_COEFF/2)
                mdfint = np.roll(mdfint, -delay)

            if self.fidVec.size != 0:
                mdfintFidel[self.fidVec>2] = 0
            else:
                mdfintFidel = mdfint

            if self.num_samples/self.fs > 90:
                xs=np.sort(mdfintFidel[self.fs:self.fs*90])
            else:
                xs= np.sort(mdfintFidel[self.fs:])
            
            if self.MAX_FORCE.size != 0:
                en_thres = self.MAX_FORCE
            else:
                if self.num_samples/self.fs>10:
                    ind_xs = ceil(98/100*xs.size)
                    en_thres = xs[ind_xs] # if more than ten seconds of ecg then 98% CI
                else:
                    ind_xs = ceil(99/100*xs.size)
                    en_thres = xs[ind_xs]

            poss_reg_bool =mdfint>(self.THRES * en_thres)

            if poss_reg_bool.size==0:
                poss_reg_bool[9] = 1

            if self.SEARCH_BACK:
                indAboveThreshold = np.where(poss_reg_bool)[0]
                RRv = np.diff(self.time[indAboveThreshold])
                medRRv = np.median(RRv[np.where(RRv>0.01)[0]])
                indMissedBeat = np.where(RRv>(1.5*medRRv))[0]
                indStart = indAboveThreshold[indMissedBeat]
                indEnd = indAboveThreshold[indMissedBeat+1]

                for i in range(0,len(indStart)):
                    poss_reg_bool[indStart[i]:indEnd[i]] = mdfint[indStart[i]:indEnd[i]] > (0.5*self.THRES*en_thres)

            poss_reg = np.zeros([poss_reg_bool.size])
            for i in range(0,poss_reg_bool.size):
                poss_reg[i] = int(poss_reg_bool[i]==True)

            left = np.where(np.diff(np.append([0],poss_reg)) == 1)[0]
            right = np.where(np.diff(np.append(poss_reg,[0])) == -1)[0]

            if self.signForce.size != 0:
                sign = self.signForce[0]
            else:
                nb_s = len(left<30*self.fs)
                loc = np.zeros([nb_s])
                for j in range(0,nb_s):
                    loc[j] = int(np.argmax(abs(self.filt_ecg[left[j]:right[j]+1])) -1 +left[j])
                    #loc[j] = int(loc[j] -1 + left[j])

                tot=0
                for j in range(0,len(loc)):
                    tot = tot + self.ecg[int(loc[j])]
                sign = tot/(len(loc))
            
            #Loop through all possibilities
            compt = 0
            NB_PEAKS = np.size(left)
            maxval = np.array([])
            maxloc = np.array([])
            for i in range(0,NB_PEAKS):
                if sign > 0:
                    maxloc = np.append(maxloc, np.argmax(self.ecg[left[i]:right[i]+1]))
                    maxval = np.append(maxval, self.ecg[int(maxloc[compt])])
                else:
                    maxloc = np.append(maxloc, np.argmin(self.ecg[left[i]:right[i]+1]))
                    maxval = np.append(maxval, self.ecg[int(maxloc[compt])])
                maxloc[compt] = maxloc[compt] + left[i] # add offset of present location

            # Refractory period --> improves results
                if compt > 0:
                    if (maxloc[compt] - maxloc[compt-1] < (self.fs * self.refPeriod)) and (abs(maxval[compt]) < abs(maxval[compt-1])):
                        maxloc = np.delete(maxloc,compt)
                        maxval = np.delete(maxval,compt)
                    elif (maxloc[compt] - maxloc[compt-1] < (self.fs * self.refPeriod)) and (abs(maxval[compt]) >= abs(maxval[compt-1])):
                        maxloc = np.delete(maxloc,compt-1)
                        maxval = np.delete(maxval,compt-1)
                    else:
                        compt=compt+1
                else:
                    compt=compt+1

            self.qrs_pos = maxloc # datapoints qrs positions
            self.R_t = np.zeros([maxloc.size]) #timestamps QRS positions
            for i in range(0,maxloc.size):
                self.R_t[i] = self.time[int(maxloc[i])] 
            self.R_amp = maxval # amplitude at QRS positions
            self.hrv = 60*np.ones([len(self.R_t)-1]) / (np.diff(self.R_t))
            self.hrv = np.append([0],self.hrv)

            return self.qrs_pos, sign, en_thres

        def plotDebug(self):
            fig, ax = plt.subplots(3, 1, figsize = (30,6), sharex=True, num = 0)
            ax[0].plot(self.time, self.ecg, 'b', zorder = 0, label="raw ecg")
            #ax[0].plot(self.time, self.filt_ecg, 'r', zorder=5, label="filtered ecg")
            #ax[0].legend()

            ax[1].plot(self.time,self.filt_ecg)
            ax[1].plot(self.R_t,self.R_amp, '+k')

            ax[2].plot(self.R_t[:len(self.hrv)], self.hrv, 'r+')
            ax[2].set_ylim([0,200])
            plt.show()
            plt.close()

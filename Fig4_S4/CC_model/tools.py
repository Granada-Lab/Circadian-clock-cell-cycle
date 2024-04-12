import numpy as np
from matplotlib import pyplot as plt
#  plt.style.use('science')
from scipy import stats

def master_sort(CL):
    N = CL.shape[0]
    CL_max = CL.max(axis = 1)
    pro_set = np.array(list(set(CL_max))).astype(int)
    key = []
    for pro in pro_set:
        idx = np.where(CL_max == pro)[0]
        if pro == 0:
            #  FP_index = idx.copy()
            #  key = list(idx)
            key.append(idx)
            continue

        FP_index = []
        for index in idx:
            FP = np.nonzero(CL[index,:])[0][0]
            FP_index.append(FP)
        FP_index = np.array(FP_index)
        FP_key = np.argsort(FP_index)
        FP_index = FP_index[FP_key]
        idx = idx[FP_key]
        key.append(idx)
    key = np.hstack(key)
    return key

def extract_cell_division(CL):
    N = CL.shape[0]
    CL_max = CL.max(axis = 1)
    pro_set = np.array(list(set(CL_max))).astype(int)
    lib = {}
    lib_time = {}
    total_peaks = []
    total_times = []
    for pro in pro_set:
        idx = np.where(CL_max == pro)[0]
        if pro < 2:
            continue
        peak_list = []
        time_list = []
        for index in idx:
            dCL = np.diff(CL[index,:])
            peaks = np.where(dCL == 1)[0]
            dPeaks = np.diff(peaks)
            peak_list.append(dPeaks)
            time_list.append(peaks[:-1])
        peak_list = np.hstack(peak_list)
        time_list = np.hstack(time_list)
        total_peaks.append(peak_list)
        total_times.append(time_list)
        lib[f"prolif_{pro}"] = peak_list
        lib_time[f"prolif_{pro}"] = time_list
    total_peaks = np.hstack(total_peaks)
    total_times = np.hstack(total_times)
    lib["prolif_total"] = total_peaks
    lib_time["prolif_total"] = total_times

    return lib, lib_time

def cells_correlation(CL):
    period_lib = {}
    times_lib = {}
    CL_max = CL.max(axis = 1)
    #  print(CL_max)
    for i in range(1,int(np.max(CL_max)),1):
        period_lib[f'div_{i}'] = []
        times_lib[f'div_{i}'] = []

    #  print(CL.shape)

    for i in range(CL.shape[0]):
        if CL_max[i] < 2:
            continue
        diff = np.diff(CL[i,:])
        idx = np.where(diff == 1)[0]
        for j in range(len(idx)-1):
            period = idx[j+1] - idx[j]
            period_lib[f'div_{j+1}'].append(period)
            times_lib[f'div_{j+1}'].append(idx[j])
    #  for i in range(1,int(np.max(CL_max)),1):
        #  if len(times_lib[f'div_{i}']) < 1:
            #  continue
        #  print(stats.pearsonr(times_lib[f'div_{i}'], period_lib[f'div_{i}'])[0])
    return times_lib, period_lib








def prolif_matrix(CL):
    dCL = np.diff(CL, axis = 1)
    sumCL = np.sum(dCL, axis = 0)
    return np.cumsum(sumCL)

def run_mean(x, N):
    return np.convolve(x, np.ones(N)/N, mode='valid')




import numpy as np
from matplotlib import pyplot as plt
plt.style.use('K_PAPER')
from tools import *
from scipy.signal import find_peaks
import tqdm
from appstatpy.ExternalFunctions import *
from scipy import stats, special
from iminuit import Minuit
from poincare_oscillator import coupling
from multiprocessing import Pool
#  output = '/home/nordentoft/Nextcloud/Manuscripts/Circadian-CellCycle/Malthe_figures/'
'''
author @malthenielse
malthe.nielsen@nbi.ku.dk
'''

def exponnorm(x, mu, sigma, lam, N):
    erf = 1 - special.erf((mu + lam*sigma**2 - x)/(np.sqrt(2)*sigma))
    return N * lam/2*np.exp(lam/2*(2*mu + lam*sigma**2 - 2*x))*erf

def log2(x, xhalf, beta):
    return 1/(1 + np.exp(beta*(x - xhalf)))

def run(circ, period):
    circ = np.roll(circ, np.random.randint(30,200,1)[0])
    RT = 0
    circ = (1/2 + 1/2*circ)**1.6
    #  circ /= np.max(circ)
    circ = np.nan_to_num(circ)

    Vol = 2

    G1 = 3
    G2 = 0 
    S = 0
    Div = 0
    AM = 0
    #  RT = 0
    Sav = np.zeros((100000, 6))

    k1 = 50
    k2 = 20
    k3 = 200
    k4 = 9
    k5 = 0.7
    k6 = 0.5
    k6 = 0.9
    #  k6 = 0.01
    k7 = 1
    k8 = 5e-3

    alpha = 0
    #  alpha = 0

    GMAX = period*3.8 + stats.expon.rvs(0, 1/.005,1)[0]
    GMAX = int(GMAX)



    rSav = np.zeros((4,280))

    for i in range(100000):

        ptime = int((RT/.1))

        dG1p = k1*(1-S)*Vol
        dG1n = k2 * G1 / (G1 + k3)*Vol
        dS = k5*max(G1-GMAX*Vol,0)*(1-S)*circ[ptime]

        dG2p = k4*S*Vol 
        dDiv = k6 * log2(G2, 150*Vol, -1)*Vol#*circ[ptime]
        dAM = (1-S)*k7*Vol
        dA  = 0*k8*log2(AM, 10*Vol, -1)*Vol

        rates = np.array([dG1p, dG1n, dG2p, dDiv, dS, dAM, dA])
        trates = np.sum(rates)
        DT = 1/trates * np.log(1/np.random.uniform(0,1,1)[0])
        rates /= trates

        ran = np.random.uniform(0,1,1)[0]
        tsum = 0
        idx = 0
        for j in range(7):
            tsum += rates[j]
            if ran < tsum:
                idx = j
                break

        if idx == 0:
            G1 += 1

        elif idx == 1:
            G1 -= 1

        elif idx == 2:
            G2 += 1

        elif idx == 3:
            Div += 1
            G2 = 0
            S = 0
            rSav[0,int(RT*2):] = Div

        elif idx == 4:
            S = 1
            AM = 1

        elif idx == 5:
            AM += 1

        elif idx == 6:
            Sav[i:, :] = np.nan
            rSav[-1,int(RT*2)] = 1

            break

        Sav[i,0] = G1
        Sav[i,1] = G2
        Sav[i,2] = circ[ptime]
        Sav[i,3] = RT
        Sav[i,4] = dA
        Sav[i,5] = AM
        rSav[1,int(RT*2)] = G1
        rSav[2,int(RT*2)] = G2


        RT += DT

        if RT > 140:
            Sav[i:, :] = np.nan
            break

    return Sav, rSav

def runner(idx):
    kinds = ['none', 'medium', 'high']
    kind = kinds[idx]
    print(kind)
    data = np.genfromtxt(f'output/circ_{kind}.csv', delimiter = ',')
    period_cc = np.genfromtxt(f'output/circ_dist_{kind}.csv', delimiter = ',')
    data = data[:,:-1]
    #  data = np.roll(data, -100, axis = 1)

    cells = np.zeros((2, 800,280))
    for i in tqdm.tqdm(range(0,800,1)):
        Sav, rSav =  run(data[i,:], period_cc[i])
        cells[:, i,:] = rSav[[0,-1],:]
    np.save(f'final_cells/cells_{kind}_S', cells)


if __name__ == '__main__':
    index = np.arange(0,3,1)
    pool = Pool(3)
    res = pool.map(runner, index)
    pool.close()
    pool.join()





#  #  coupling(0.000001, 'none')
#  data = np.genfromtxt('nica_output/circ_none.csv', delimiter = ',')
#  period_cc = np.genfromtxt('nica_output/circ_dist_none.csv', delimiter = ',')
#  data = data[:,:-1]
#  data = np.roll(data, -100, axis = 1)
#
#  cells = np.zeros((2, 800,280))
#  for i in tqdm.tqdm(range(0,800,1)):
#      Sav, rSav =  run(data[i,:], period_cc[i])
#      cells[:, i,:] = rSav[[0,-1],:]
#  np.save('final_cells/cells_none_G1', cells)
#
#  #  coupling(0.005, 'high_alt')
#  data = np.genfromtxt('nica_output/circ_medium.csv', delimiter = ',')
#  period_cc = np.genfromtxt('nica_output/circ_dist_medium.csv', delimiter = ',')
#  data = data[:,:-1]
#  data = np.roll(data, -100, axis = 1)
#  cells = np.zeros((2, 800,280))
#  for i in tqdm.tqdm(range(800)):
#      Sav, rSav =  run(data[i,:], period_cc[i])
#      cells[:, i,:] = rSav[[0,-1],:]
#  np.save('final_cells/cells_medium_G1', cells)
#
#  #  coupling(0.005, 'high_alt')
#  data = np.genfromtxt('nica_output/circ_high.csv', delimiter = ',')
#  period_cc = np.genfromtxt('nica_output/circ_dist_high.csv', delimiter = ',')
#  data = data[:,:-1]
#  data = np.roll(data, -100, axis = 1)
#  cells = np.zeros((2, 800,280))
#  for i in tqdm.tqdm(range(800)):
#      Sav, rSav =  run(data[i,:], period_cc[i])
#      cells[:, i,:] = rSav[[0,-1],:]
#  np.save('final_cells/cells_high_G1', cells)

    





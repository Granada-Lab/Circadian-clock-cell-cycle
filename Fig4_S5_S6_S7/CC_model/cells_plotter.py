import numpy as np
from matplotlib import pyplot as plt
plt.style.use('K_PAPER')
from tools import *
from appstatpy.ExternalFunctions import *
from scipy import stats, special
from iminuit import Minuit
from scipy.optimize import curve_fit
#  output = '/home/nordentoft/Nextcloud/Manuscripts/Circadian-CellCycle/Malthe_figures/'
from matplotlib.colors import LinearSegmentedColormap

def create_custom_colormap(num_colors): #input: maximum number divisions
    all_colors = [
        (135/255, 206/255, 250/255),  # skyblue
        (255/255, 228/255, 181/255),  # moccasin
        (216/255, 191/255, 216/255),  # thistle
        (240/255, 128/255, 128/255),  # lightcoral
        (95/255, 158/255, 160/255),   # cadetblue
        (255/255, 215/255, 0/255),    # gold
    ]

    if num_colors not in range(3, 7):
        raise ValueError("Number of colors must be between 3 and 6 inclusive")

    colors = all_colors[:num_colors]

    cmap_name = 'custom_colormap'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

ccmap = create_custom_colormap(6)

def exponnorm(x, mu, sigma, lam, N):
    erf = 1 - special.erf((mu + lam*sigma**2 - x)/(np.sqrt(2)*sigma))
    return N * lam/2*np.exp(lam/2*(2*mu + lam*sigma**2 - 2*x))*erf

def lin_func(x, a, b):
    return a*x + b

#  kind = 'S'
kind = 'G1'
#  kind = 'G2'
#  kind = 'M'

ch = np.load(f'final_cells/cells_high_{kind}.npy')
cn = np.load(f'final_cells/cells_none_{kind}.npy')
cm = np.load(f'final_cells/cells_medium_{kind}.npy')

dch = ch[1,:,:]
ch  = ch[0,:,:]
dcn = cn[1,:,:]
cn  = cn[0,:,:]
dcm = cm[1,:,:]
cm  = cm[0,:,:]

def div_matrix(mat, string):
    dmat = np.diff(mat, axis = 1)
    dmat = dmat[:,40:]
    np.savetxt('dmats/'+string+'.csv', dmat, delimiter = ',')
#
div_matrix(ch, 'high' + '_' + kind )
div_matrix(cn, 'none' + '_' + kind)
div_matrix(cm, 'midle' + '_' + kind)


#  exit()

#  ch = np.load('./nica_output/pmat_high.npy')
#  cn = np.load('./nica_output/pmat_none.npy')
#  ch = np.cumsum(ch, axis = 1)
#  cn = np.cumsum(cn, axis = 1)


#  dch = np.cumsum(np.sum(dch, axis =0))
#  dcn = np.cumsum(np.sum(dcn, axis =0))


#  time = np.arange(20.5, ch.shape[1]//2 - 40, .5)


#  def sin_lin(time, alpha, beta, A0, T, phase):
#      return alpha*time + beta + A0*np.sin(2*np.pi/T*time + phase)
#
#  pm = prolif_matrix(ch)/400
#  pm = pm[55:]
#  time = np.arange(0, len(pm)/2, .5)
#  fig, ax = plt.subplots(1,2, figsize = (14,7))
#  ax[0].plot(time, pm, color = 'black')
#  p0, _ = curve_fit(sin_lin, time, pm, p0 = [1/24, 0, 10, 20, 0])
#  ax[0].plot(time, sin_lin(time, *p0), color = 'black', ls = '--')
#  print(p0, 1/p0[0])
#
#  ax[1].plot(time, pm - (p0[0]*time + p0[1]))
#  ax[1].plot(time, p0[2]*np.sin(2*np.pi/p0[3]*time + p0[4]))
#
#  plt.show()
#  exit()

#  fig, ax = plt.subplots(1,1, figsize = (7,7))




idx = master_sort(ch)
ch = ch[idx, :]

idx = master_sort(cn)
cn = cn[idx, :]

idx = master_sort(cm)
cm = cm[idx, :]

fig, ax = plt.subplots(1,3, figsize = (7,10), sharex = True, sharey = True)
ax[0].imshow(ch, aspect = 'auto', extent = [-20,120, 0, 400], cmap = ccmap, vmin = 0, vmax = 5)
ax[1].imshow(cm, aspect = 'auto',extent = [-20,120, 0, 400] , cmap = ccmap, vmin = 0, vmax = 5)
ax[2].imshow(cn, aspect = 'auto',extent = [-20,120, 0, 400] , cmap = ccmap, vmin = 0, vmax = 5)
ax[0].set_xlim(0,120)
ax[0].set(title = 'Untreated', xlabel = 'Time [hr]', ylabel = 'Cell #')
ax[1].set(title = '5uM', xlabel = 'Time [hr]')
ax[2].set(title = '10uM', xlabel = 'Time [hr]')
fig.savefig(output + 'prolif_matrix_simulation.svg', dpi = 500)

#  plt.show()

time = np.arange(.5, ch.shape[1]/2, .5)
#  print(ch.shape)


fig, ax = plt.subplots(1,1, figsize = (7,5))
pch = prolif_matrix(ch)
pcm = prolif_matrix(cm)
pcn = prolif_matrix(cn)
max_val = max([pch[-1], pcm[-1], pcn[-1]])
max_val = ch.shape[0]

grate = []

ax.plot(time[40:], pch[40:]/max_val, color = 'black', label = 'High')
p0, _ = curve_fit(lin_func, time[40:],prolif_matrix(ch)[40:]/max_val)
p00, _ = curve_fit(lin_func, time[40:],prolif_matrix(ch)[40:]/800)
grate.append(p0[0])
#  print(time[40])
#  print(p0)
#  ax.plot(time[40:], lin_func(time[40:], *p0), color = 'black', ls = '--', label = f'Slope = {np.round(1/p00[0],1)}')

ax.plot(time[40:], pcm[40:]/max_val, color = 'red', label = 'low')
p0, _ = curve_fit(lin_func, time[40:], prolif_matrix(cm)[40:]/max_val)
p00, _ = curve_fit(lin_func, time[40:], prolif_matrix(cm)[40:]/800)
grate.append(p0[0])
#  ax.plot(time[40:], lin_func(time[40:], *p0), color = 'red', ls = '--', label = f'Slope = {np.round(1/p00[0],1)}')

ax.plot(time[40:], pcn[40:]/max_val, color = 'green', label = 'None')
p0, _ = curve_fit(lin_func, time[130:],prolif_matrix(cn)[130:]/max_val)
p00, _ = curve_fit(lin_func, time[130:],prolif_matrix(cn)[130:]/800)
grate.append(p0[0])
#  print(1/p0)
#  ax.plot(time[40:], lin_func(time[40:], *p0), color = 'green', ls = '--', label = f'Slope = {np.round(1/p00[0],1)}')

ax.set_xticks(np.linspace(20,150,6), [0, 24, 48, 72, 96, 120])

ax.legend()

#  fig, ax = plt.subplots(1,1, figsize = (7,7))
#  ax[1].plot(time[14:-5], run_mean(np.diff(run_mean(prolif_matrix(ch),10)), 10), color = 'black')
#  ax[1].plot(time[14:-5], run_mean(np.diff(run_mean(prolif_matrix(cn),10)), 10), color = 'green')
#
#  ax[1].set(xlabel = 'Time [hr]', ylabel = 'Signal[AU]', title = 'Detrended growth')
ax.set(xlabel = 'Time [hr]', ylabel = 'Total fate', title = 'Growth')
fig.savefig(output + 'Growth_simulation.svg', dpi = 500)



fit = {'ExpNorm' : ['$\mu$', '$\sigma_{\mu}$']}
fig, ax = plt.subplots(1,1, figsize = (7,7))
bins = np.arange(10, 100,1)

err_bin = []

res = np.zeros((3,3))

lib, _ = extract_cell_division(ch)
#  lib['prolif_total'] = np.delete(np.where(lib['prolif_total'] > 70)[0],lib['prolif_total'])
lib['prolif_total'] = lib['prolif_total'][np.where(lib['prolif_total'] < 90)[0]]

#  print(lib['prolif_total'])
#  ax.hist(lib['prolif_total']/2, bins = bins, histtype = 'step', lw = 2, density = True, color = 'black')
ax.hist(np.array(lib['prolif_total'])/2, histtype = 'step', color = 'black', bins = np.arange(10,60,.5), density = True)
freq, bins = np.histogram(lib['prolif_total'], bins = np.arange(10,60,.5), density = False)
freq_err = np.sqrt(freq)
bin_center = (bins[:-1] + bins[1:])/4
bw = bin_center[1] - bin_center[0]
chi_im = Chi2Regression(exponnorm, bin_center[freq > 0], freq[freq > 0], freq_err[freq > 0])
chi_im = Minuit(chi_im, mu = 20, sigma = 2.1, lam = 0.1, N = np.sum(freq)*bw)
chi_im.errordef = 1
chi_im.migrad()
ullh_x = np.linspace(0,100,1000)
print(chi_im.values[:])
ax.plot(ullh_x, exponnorm(ullh_x, *chi_im.values[:-1], N = 1), color = 'black')
mu_true = chi_im.values[0] + 1/chi_im.values[2]
chi2 = stats.chi2.sf(chi_im.fval, len(freq[freq >0]))
err = np.sqrt(chi_im.errors[0]**2 + (1/chi_im.values[2]**2)*chi_im.errors[2]**2)
print(mu_true, err, chi2)
res[0,:] = chi_im.values[:3]
err_bin.append(err)

fit['Untreated'] = [mu_true, err]

lib, _ = extract_cell_division(cm)
lib['prolif_total'] = lib['prolif_total'][np.where(lib['prolif_total'] < 90)[0]]
ax.hist(np.array(lib['prolif_total'])/2, histtype = 'step', color = 'red', bins = np.arange(10,60,.5), density = True)
freq, bins = np.histogram(lib['prolif_total'], bins = np.arange(10,60,.5), density = False)
freq_err = np.sqrt(freq)
bin_center = (bins[:-1] + bins[1:])/4
bw = bin_center[1] - bin_center[0]
chi_im = Chi2Regression(exponnorm, bin_center[freq > 0], freq[freq > 0], freq_err[freq > 0])
chi_im = Minuit(chi_im, mu = 20, sigma = 2.1, lam = 0.1, N = np.sum(freq)*bw)
chi_im.errordef = 1
chi_im.migrad()
ullh_x = np.linspace(0,100,1000)
print(chi_im.values[:])
ax.plot(ullh_x, exponnorm(ullh_x, *chi_im.values[:-1], N = 1), color = 'red')
mu_true = chi_im.values[0] + 1/chi_im.values[2]
chi2 = stats.chi2.sf(chi_im.fval, len(freq[freq >0]))
err = np.sqrt(chi_im.errors[0]**2 + (1/chi_im.values[2]**2)*chi_im.errors[2]**2)
print(mu_true, err, chi2)
res[1,:] = chi_im.values[:3]
fit['5uM'] = [mu_true, err]
err_bin.append(err)

lib, _ = extract_cell_division(cn)
lib['prolif_total'] = lib['prolif_total'][np.where(lib['prolif_total'] < 90)[0]]
ax.hist(np.array(lib['prolif_total'])/2, histtype = 'step', color = 'green', bins = np.arange(10,60,.5), density = True)
freq, bins = np.histogram(lib['prolif_total'], bins = np.arange(10,60,.5), density = False)
freq_err = np.sqrt(freq)
bin_center = (bins[:-1] + bins[1:])/4
bw = bin_center[1] - bin_center[0]
chi_im = Chi2Regression(exponnorm, bin_center[freq > 0], freq[freq > 0], freq_err[freq > 0])
chi_im = Minuit(chi_im, mu = 18, sigma = 2.3, lam = 0.1, N = np.sum(freq)*bw)
chi_im.errordef = 1
chi_im.migrad()
ullh_x = np.linspace(0,100,1000)
print(chi_im.values[:])
ax.plot(ullh_x, exponnorm(ullh_x, *chi_im.values[:-1], N = 1), color = 'green')
mu_true = chi_im.values[0] + 1/chi_im.values[2]
chi2 = stats.chi2.sf(chi_im.fval, len(freq[freq >0]))
err = np.sqrt(chi_im.errors[0]**2 + (1/chi_im.values[2]**2)*chi_im.errors[2]**2)
print(mu_true, err, chi2)
err_bin.append(err)

res[2,:] = chi_im.values[:3]
fit['10uM'] = [mu_true, err]
#  ax.plot(ullh_x, exponnorm(ullh_x, 24, 5, .1, N = 1), color = 'red')

grate = np.array([grate])

res = np.concatenate((res, grate), axis = 0)
np.save(f'res_{kind}.npy', res)

text = nice_string_output(fit, extra_spacing=0, decimals=3)
add_text_to_ax(.5, .5, text, ax, fontsize = 11)
ax.set(xlabel = 'Period [hr]', ylabel = 'p', title = 'Intermiotic time', xlim = (15, 55))
#  fig.savefig(output + 'IM_historgram_simulation.svg', dpi = 500)

#  time = np.arange(0, len(avg)/2, .5)
#  axx[0].plot(time, avg/max_val, color = colors[i], alpha = 1)
#  axx[0].plot(time, avg/351, color = colors[i], alpha = 1)

#  p0, _ = curve_fit(lin_func, time[:120], avg[:120]/max_val, p0 = [1, 1])
#  p0, _ = curve_fit(lin_func, time[:90], avg[:90]/351, p0 = [1, 1])


plt.show()


print(err_bin)




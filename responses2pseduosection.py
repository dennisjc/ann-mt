import numpy as np
from keras.models import load_model
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import random
import matplotlib.colors as colors
from matplotlib import ticker
from matplotlib_scalebar.scalebar import ScaleBar
from llres_yes import scale_params
from create_test import get_model
import subprocess

def scale_loc(x):
    minx = -0.234375E+05
    maxx = 0.234375E+05
    return (x - minx)/ (maxx - minx)

# blob = [0.8918586867152596,
#  0.6973829053990582,
#  0.558853236019389,
#  0.6253231058720568,
#  0.4985224386549934,
#  0.49283782183812214,
#  0.1189641343484739,
#  0.08293705606611435,
#  0.20183854624660633,
#  0.10755224411830816]

# blob = [0.9, 0.0, 0.01, 0.505, 0.5, 0.5, 0.125, 0.125, 0.2, 0.125]
#blob = [0.9, 0.2, 0.2, 0.505, 0.5, 0.5, 0.1, 0.125, 0.225, 0.125]
# blob = [0.82, 0.206, 1, 0.51, 0.5, 0.5, 0.065, 0.069, 0.09, 0.069]
blob = [0.9, 0.2, 0.2, 0.505, 0.5, 0.5, 0.1, 0.125, 0.225, 0.125]
blob = [0.9, 0.0, 0.01, 0.505, 0.5, 0.5, 0.125, 0.125, 0.2, 0.125]
blob = scale_params(blob).tolist()
# blobs = (blob).tolist()
get_model(blob, 1)
p1 = subprocess.Popen(['./wsinv3dmt'])
p1.wait()
#blobs = [ 0.85961756,  0.53216788,  0.01000001,  0.99      ,  0.50065602,
#        0.50196512,  0.18937985,  0.19521973,  0.23694558,  0.2       ]
# blobs = [float(i) for i in blobs.split(' ') if float(i) != 0]
#$ blobs = [i for i in blobs if i != 0]
#blobs.pop(-1)
#blobs.pop(-1)
#blobs.pop(-3)
model = load_model('oneblob.model')
mins = np.load('oneblob.min')
maxs = np.load('oneblob.max')

with open('./fwd_resp.02') as f:
    data = f.readlines()

header = [float(i) for i in data.pop(0).split(' ') if i]
n_sites, freqs, components = header[0], header[1], header[2]
data.pop(0)
x_locs = []

while True:
    line = data.pop(0)
    if 'Station' in line:
        break
    x_locs.extend([float(i) for i in line.split(' ') if i])

y_locs = []
response = {}
while True:
    line = data.pop(0)
    if 'DATA' in line:
        period = float(line.split(' ')[-1])
        response[period] = []
        break
    y_locs.extend([float(i) for i in line.split(' ') if i])

while True:
    line = data.pop(0)
    if 'Iteration' in line:
        break
    if 'DATA' in line:
        period = float(line.split(' ')[-1])
        response[period] = []
        continue
    tensor = [float(i.strip()) for i in line.split(' ') if i.strip()][2:6]
    app_xy = np.linalg.norm(np.array(tensor[0:2])*796)**2*period*0.2
    phs_xy = np.degrees(np.arctan2(tensor[0],-tensor[1]))
    app_yx = np.linalg.norm(np.array(tensor[2:4])*796)**2*period*0.2
    phs_yx = np.degrees(np.arctan2(-tensor[2],tensor[3]))
    response[period].append([app_xy, phs_xy, app_yx, phs_yx])

ins = []
x_set = sorted(list(set(x_locs)))
y_set = sorted(list(set(y_locs)))
periods = sorted(response.keys())
data_matrix = np.zeros([len(x_set), len(y_set), len(response.keys()), 4])
neural_matrix = np.zeros(data_matrix.shape)
for idx, (i, j) in enumerate(zip(x_locs, y_locs)):
    x = np.argmax(np.array(x_set) == i)
    y = np.argmax(np.array(y_set) == j)
    neural_in = np.array([scale_loc(x_set[x]), scale_loc(y_set[y])] + blob).reshape(1, 12)
    ins.append(scale_loc(x))
    neural_out = (model.predict(neural_in)[0]*(maxs-mins)+mins).reshape(2, 5, 2)
    neural_out[1] = 10**neural_out[1]
    neural_out[0] = 90 - neural_out[0]
    for jdx, k in enumerate(sorted(response.keys())):
        data_matrix[x, y, jdx] = response[k][idx]
        neural_matrix[x, y, jdx] = np.array([neural_out[1, jdx, 0], neural_out[0, jdx, 0], neural_out[1, jdx, 1], neural_out[0, jdx, 1]])


f, axes = plt.subplots(3, 4, sharey=True, figsize=(6, 5.1))
for idx, row in zip([0, 2, 4], axes):
    for jdx, (nax, dax) in enumerate(zip(row[:4][::2], row[:4][1::2])):
        t = np.append(neural_matrix[:, :, :, jdx], data_matrix[:, :, :, jdx])
        if jdx == 0:
            # continue
            im1 = nax.pcolormesh(x_set, y_set, np.log10(neural_matrix[:, :, idx, jdx]), cmap='viridis_r', edgecolor='k', linewidth=0.005, antialiased=False)
            im2 = dax.pcolormesh(x_set, y_set, np.log10(data_matrix[:, :, idx, jdx]), cmap='viridis_r', edgecolor='k', linewidth=0.005, antialiased=False)
            im1.set_clim([np.log10(t.min()), np.log10(t.max())])
            im2.set_clim([np.log10(t.min()), np.log10(t.max())])
            cbar_ax = f.add_axes([0.125, 0.19, 0.5-0.1275, 0.0333333])
            cb = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
            cb.set_label(r'$log_{10}[\rho_{a,xy}]$ $(\Omega m)$')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            cb.update_ticks()
            # cb.ax.tick_params(axis='x',direction='in',labeltop='on')
        else:
            im1 = nax.pcolormesh(x_set, y_set, neural_matrix[:, :, idx, jdx], cmap='magma', edgecolor='k', linewidth=0.005, antialiased=False)
            im2 = dax.pcolormesh(x_set, y_set, data_matrix[:, :, idx, jdx], cmap='magma', edgecolor='k', linewidth=0.005, antialiased=False)
            im1.set_clim([t.min(), t.max()])
            im2.set_clim([t.min(), t.max()])
            cbar_ax = f.add_axes([0.53, 0.19, 0.5-0.1275, 0.0333333])
            cb = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
            cb.set_label(r'$\phi_{xy}$')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cb.locator = tick_locator
            cb.update_ticks()
            # cb.ax.tick_params(axis='x',direction='in',labeltop='on')
        nax.set_xlim([x_set[0], x_set[-1]])
        dax.set_xlim([x_set[0], x_set[-1]])
        nax.set_ylim([y_set[0], y_set[-1]])
        dax.set_ylim([y_set[0], y_set[-1]])
        if jdx != 0:
            nax.tick_params(labelleft='off')
        else:
            nax.set_ylabel('T = {} s'.format(periods[idx]).replace('0.56', '0.5'))
            nax.tick_params(labelleft='off')
            nax.set_xticklabels(['-2 km', '', '0 km', '', '2 km'])
        if idx != 0:
            nax.tick_params(labelbottom='off')
            dax.tick_params(labelbottom='off')
        else:
            bot = r'{} ${}_{{xy}}$'
            # nax.set_xlabel(bot.format('ANN', '\phi' if jdx == 1 else r'\rho'))
            nax.set_xlabel('ANN')
            nax.tick_params(labelbottom='off')
            nax.xaxis.set_label_position('top')
            # dax.set_xlabel(bot.format('WSINV3D', '\phi' if jdx == 1 else r'\rho'))
            dax.set_xlabel('WSINV3D')
            dax.tick_params(labelbottom='off')
            dax.xaxis.set_label_position('top')
        dax.tick_params(labelleft='off')
    # ax = row[4]
    # im1 = ax.pcolormesh(x_set, y_set, neural_matrix[:, :, idx, 0] - data_matrix[:, :, idx, 0], cmap='RdBu', edgecolor='k', linewidth=0.005, antialiased=False)
    # maxi = np.max(np.abs(neural_matrix[:, :, idx, 0] - data_matrix[:, :, idx, 0]))
    # im1.set_clim([-maxi, maxi])
    # cbar_ax = f.add_axes([0.53, 0.19, 0.5-0.1275, 0.0333333])
    # cb = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
    # cb.set_label('wah')
    # tick_locator = ticker.MaxNLocator(nbins=5)
    # cb.locator = tick_locator
    # cb.update_ticks()
    # ax.set_ylim(y_set[0], y_set[-1])
    # ax.set_xlim(x_set[0], x_set[-1])
    # ax = row[5]
    # im1 = ax.pcolormesh(x_set, y_set, neural_matrix[:, :, idx, 1] - data_matrix[:, :, idx, 1], cmap='RdBu', edgecolor='k', linewidth=0.005, antialiased=False)
    # maxi = np.max(np.abs(neural_matrix[:, :, idx, 1] - data_matrix[:, :, idx, 1]))
    # im1.set_clim([-maxi, maxi])
    # cbar_ax = f.add_axes([0.53, 0.19, 0.5-0.1275, 0.0333333])
    # cb = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
    # cb.set_label('wah')
    # tick_locator = ticker.MaxNLocator(nbins=5)
    # cb.locator = tick_locator
    # cb.update_ticks()
    # ax.set_ylim(y_set[0], y_set[-1])
    # ax.set_xlim(x_set[0], x_set[-1])

sb_ax = f.add_axes([0.125, 0.02, 1.-0.125*2+0.025+0.02, 0.05])
sb_ax.axes.get_xaxis().set_visible(False)
sb_ax.axes.get_yaxis().set_visible(False)
sb_ax.axis('off')
scalebar = ScaleBar(175000.)
plt.gca().add_artist(scalebar)
f.subplots_adjust(bottom=0.25)

plt.savefig('./forward_rev.png', dpi=1000)
plt.savefig('./forward_rev.pdf')

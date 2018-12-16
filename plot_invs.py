import numpy as np
from keras.models import load_model
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import random
import matplotlib.colors as colors
from matplotlib import ticker
from ws2vtk import main as ws2vtk
from collections import OrderedDict
from matplotlib_scalebar.scalebar import ScaleBar

ann = ws2vtk('ann.model', '../oneblob.resp')
# inv = ws2vtk('targ2.model', '../oneblob.resp')
inv = ws2vtk('blobws.model', '../oneblob.resp')
truth = ws2vtk('truth.model', '../oneblob.resp')

models = OrderedDict()
models['True model'] = truth
models['ANN inversion'] = ann
models['WSINV3D inversion'] = inv
labels = {4: '1906-2431 m', 8: '4059-4616 m', 12: '6317-6898 m', 16: '8646-9237 m', 18: '9831-10430 m'}

x_set = np.linspace(-25, 25, 20)
y_set = np.linspace(-25, 25, 20)

all = np.array([truth, ann, inv])

f, axes = plt.subplots(4, 3, sharey=True, figsize=(6, 8.5))

for idx, row in zip([4, 8, 12, 16], axes):
    for jdx, (ax, model) in enumerate(zip(row, models)):
        im1 = ax.pcolormesh(x_set, y_set, np.log10(models[model][:, :, idx]), cmap='inferno_r', edgecolor='k', linewidth=0.005, antialiased=False)
        im1.set_clim([np.log10(all.min()), np.log10(all.max())])
        # cb.ax.tick_params(axis='x',direction='in',labeltop='on')
        ax.set_xlim([x_set[0], x_set[-1]])
        ax.set_ylim([y_set[0], y_set[-1]])
        print ax.get_xticks().tolist()
        # print labels
        # ax.set_xticklabels(labels)
        # plt.setp(ax)
        # plt.tight_layout()
        if jdx != 0:
            ax.tick_params(labelleft='off')
        else:
            ax.set_ylabel(labels[idx])
            ax.tick_params(labelleft='off')
        if idx != 4:
            ax.tick_params(labelbottom='off')
        else:
            bot = r'{} ${}_{{xy}}$'
            ax.set_xlabel(model)
            ax.tick_params(labelbottom='off')
            ax.xaxis.set_label_position('top')
        ax.tick_params(labelleft='off')
        # ax.grid(True, which='minor', axis='both', linestyle='-', color='k')
cbar_ax = f.add_axes([0.125, 0.145 - 0.0033333333, 1.-0.125*2+0.025, 0.0333333333333])
cb = f.colorbar(im1, cax=cbar_ax, orientation='horizontal')
cb.set_label(r'$log_{10}[\rho]$ $(\Omega m)$')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
sb_ax = f.add_axes([0.125, 0.025, 1.-0.125*2+0.025, 0.05])
sb_ax.axes.get_xaxis().set_visible(False)
sb_ax.axes.get_yaxis().set_visible(False)
sb_ax.axis('off')
scalebar = ScaleBar(175000.)
plt.gca().add_artist(scalebar)

f.subplots_adjust(bottom=0.2)
plt.savefig('inversions_new.pdf')
plt.savefig('inversions_new.png', dpi=1000)

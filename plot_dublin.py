import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import random
import matplotlib.colors as colors
from matplotlib import ticker
from ws2vtk import main as ws2vtk
from collections import OrderedDict
from matplotlib_scalebar.scalebar import ScaleBar

ann = ws2vtk('../dublin/dublin_ann.model', '../oneblob.resp')
inv = ws2vtk('../dublin/dublin_ws.model', '../oneblob.resp')
truth = np.load('../dublin/dublin_truth.npy')

models = OrderedDict()
models['True model'] = truth
models['ANN blob inversion'] = ann
models['FD blob inversion'] = inv
labels = ['Slice through x = 0', 'Slice through y = 0',
          'Slice through z = 30 km']

x_set = np.linspace(-25, 25, 20)
y_set = np.linspace(-25, 25, 20)

all = np.array([truth, ann, inv])

f, axes = plt.subplots(3, 3, sharey=True, figsize=(8.5, 8.5))

for idx, row in zip(labels, axes):
    for jdx, (ax, model) in enumerate(zip(row, models)):
        im1 = ax.pcolormesh(x_set, y_set, np.log10(models[model][:, :, idx]),
                            cmap='inferno_r', edgecolor='k', linewidth=0.005,
                            antialiased=False)
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
cbar_ax = f.add_axes([0.125, 0.145 - 0.0033333333, 1.-0.125*2+0.025,
                      0.0333333333333])
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
plt.savefig('dublin_inversions.pdf')
plt.savefig('dublin_inversion.png', dpi=1000)

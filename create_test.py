import sys
import random
import numpy as np
import os
import subprocess
import uuid
import shutil
import time
import pickle
from ws2vtk import main as ws2vtk
from llres_maybe import scale_params
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

## Dennis Conway

#
# with open('wae.pickle') as f:
#     b = pickle.load(f)


blob = [1, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
# background, intensity, strength, attenuation, xpos, ypos, xrad, yrad, theta, zpos, zrad, xrot, yrot
def pad_model(b, no_blobs):
    no_blobs = len(b[1:])/9
    b = list(b)
    # b[6] = 0
    inserts = []
    for i in range(no_blobs):
        for j in [7, 10, 11]:
            inserts.append(i*9+1+j)
    for i in inserts[::-1]:
        b.insert(i, 0)
    return b


def get_model(b, no_blobs):
    b = pad_model(b, no_blobs)
    c=' '.join(map(str,b))
    print(c)

    model_template = './template'
    target_model = './targ'

    callString='./saveModel '+ model_template + ' ' + c + ' > ' + target_model
    os.system(callString)
    time.sleep(2)
    # mtNS, N, E, D = ws2vtk(target_model, './oneblob.resp')

    # return mtNS, N, E, D

def test_yr_blobs(fwd, blobs):
    rms = [runw3(fwd, scale_params(blobs[0]), 1, off_diag=False)]
    for idx, b in enumerate(blobs[1:]):
        if (blobs[idx] != b).any():
            rms.append(runw3(fwd, scale_params(b), 1, off_diag=False))
        else:
            rms.append(rms[-1])
    return rms

def runw3(inv, b, no_blobs, off_diag=True):
    get_model(b, no_blobs)
    p1 = subprocess.Popen(['./wsinv3dmt'])
    p1.wait()
    #rmses = [runw3(hey, pad_model(scale_params(i[0]), 1), 1) for i in traj[:200]]
    time.sleep(12)
    # return
    #storing response and model files with unique identifier
    response='./fwd_resp.02.new'
    inv.w3 = inv.pop_data(response)
    rms = np.sqrt(np.average((inv.w3 - inv.data)**2/inv.data_errors**2))
    return rms

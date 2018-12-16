import sys
import keras
import numpy
import os, sys
import numpy as np
import shutil
import time
from glob import glob1
import cma
from create_test import runw3, pad_model
import string
import random
import os
import pickle
from multiprocessing import Pool, TimeoutError
import subprocess32 as subprocess
from dublin_inversion import scale_params, GeoFitnessNeural, evol
FNULL = open(os.devnull, 'w')


def scale_params(parameters, unscale=False):
    try:
        mins = np.array([0, 0, 0.1, 0, 0.01, 0.02, 0.4, 0.4, 0.05, 0.05, 0, 0.05])
        maxs = np.array([1, 1, 0.9, 1, 0.99, 0.99, 0.6, 0.6, 0.2, 0.2, 0.4, 0.2])
        return (parameters*(maxs-mins) + mins)
    except:
        mins = mins[2:]
        maxs = maxs[2:]
        return (parameters*(maxs-mins) + mins)


def evaluate_single(x, data, errors, xy):
    # try:
    os.chdir('/network/cogme/outputs/dennis/neural_paper')
    id_ = ''.join(random.choice(string.ascii_letters) for m in xrange(10))
    b = pad_model(scale_params(x), 1)
    c=' '.join(map(str,b))
    model_template = '/network/cogme/outputs/dennis/neural_paper/template'
    target_model = '/network/cogme/outputs/dennis/neural_paper/{}/targ'.format(id_)
    os.mkdir(id_)
    shutil.copy('saveModel', id_)
    shutil.copy('startup', id_)
    shutil.copy('wsinv3dmt', id_)
    os.chdir(id_)
    callString='./saveModel '+ model_template + ' ' + c + ' > ' + target_model
    os.system(callString)
    time.sleep(2)
    # shutil.copy('startup', id_)
    # shutil.copy('wsinv3dmt', id_)
    # os.chdir(id_)
    p1 = subprocess.Popen(['./wsinv3dmt'])
    p1.wait()
    time.sleep(2)
    model = GeoFitnessNeural(read_data=False)
    model.xy = xy
    model.data, model.data_errors = data, errors
    model.w3 = model.pop_data('./fwd_resp.02')
    rms = np.sqrt(np.average((model.w3 - model.data)**2/model.data_errors**2))
    print(rms)
    os.chdir('..')
    shutil.rmtree(id_)
    # os.rmdir(id_)
    return rms
    # except:
    #    return np.nan


def evaluate_population(x, processes, data, errors, xy):
    pool = Pool(processes=processes)
    rmses = [pool.apply_async(evaluate_single,
                              (xi, data, errors, xy)) for xi in x]
    return [i.get() for i in rmses]
    # return [evaluate_single(i, data, errors, xy) for i in x]


def evol(params, sigma, popsize, maxiter, processes=4):
    # try:
    trajectory = []
    # with open('./dublin/data.p', 'rb') as f:
    #     data, data_errors = pickle.load(f)
    fwd = GeoFitnessNeural(read_data=True)
    es = cma.CMAEvolutionStrategy(params, sigma, {
        'popsize': popsize, 'maxiter': maxiter, 'bounds': [0, 1]})
    fit = [256] * popsize
    while not es.stop():
        x = es.ask(popsize)
        fit = evaluate_population(x, processes, fwd.data, fwd.data_errors,
                                  fwd.xy)
        es.tell(x, fit)
        es.disp()
        trajectory.append(es.best.get()[:2])
    return es, trajectory
    # except:
    #     return es, trajectory

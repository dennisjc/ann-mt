import sys
import keras
import numpy
import os
import sys
import numpy as np
import time
import utm
from glob import glob1
import cma
from glob import glob
from mtpy.core.edi import Edi


def scale_loc(x):
    minx = -0.234375E+05
    maxx = 0.234375E+05
    return (float(x) - minx) / (maxx - minx)


def get_locations(filename):
    header = [], [], True
    locs = {'x': [], 'y': []}
    with open(filename) as f:
        for i, line in enumerate(f):
            line = filter(None, line.strip('\n').split(" "))
            if header:
                if 'Station' in line[0]:
                    header = False
                block = locs['x']
                continue
            if 'Station' in line[0]:
                block = locs['y']
                continue
            if 'Period' in line[0]:
                break
            block.extend([scale_loc(i) for i in line])
    return zip(locs['x'], locs['y'])


class GeoFitnessNeural():

    """ A class for reading in MT data and finding the RMS error of proposed
    blob parameters """

    def __init__(self, data_fn='blk3d.data', neural_dir='./', error=None,
                 batch_size=1024 * 4, error_floor=0.05, scale=False):
        """
            reads in data file from current directory and neural network data
            it is assumed that neural net models are structured as follows:
                base_name.model -- filename for keras hdf5 model
                base_name.min   -- filename for minimum scale parameters
                base_name.max   -- filename for maximum scale parameters

            currently only programmed for a datafile which exactly matches the
            40 parameter response files which the neural networks have been
            trained on
        """
        # self.read_edis(error_floor=error_floor)
        self.neural_nets, self.neural_scale = {}, {}
        self.batch_size = batch_size
        base_names = [i[:-6] for i in glob1(neural_dir, '*.model')]
        print base_names
        for fn in base_names:
            model = keras.models.load_model(
                os.path.join(neural_dir, fn + '.model'))
            scale_min = np.load(os.path.join(neural_dir, fn + '.min'))
            scale_max = np.load(os.path.join(neural_dir, fn + '.max'))
            if model.input_shape[1] in self.neural_nets.keys():
                print('WARNING: multiple models with the same blob numbers!')
            self.neural_nets[model.input_shape[1]] = {}
            self.neural_nets[model.input_shape[1]]['model'] = model
            self.neural_nets[model.input_shape[1]]['min'] = scale_min
            self.neural_nets[model.input_shape[1]]['max'] = scale_max
        self.read_wsinv3dmt('./dublin/RP_Miensopust_wsinv3dmt.txt')

    def read_wsinv3dmt(self, filename):
        freqs = [0.56000000000000005, 5.5999999999999996, 10.0, 100, 1000]
        with open(filename) as f:
            data = f.readlines()[1:-1]
        data = [[float(x) for x in i.strip().split(' ')[::6]] for i in data]
        fdata = [i for i in data if i[2] in freqs]
        data_matrix = np.zeros([59, 5, 4])
        error_matrix = np.zeros([59, 5, 4])
        self.x_locs, self.y_locs = [], []
        for i in range(5):
            data_matrix[:, i, 2] = [np.log10(j[5]) * (1 + np.random.randn() * 0.05)
                                    for j in fdata[i::5]]
            error_matrix[:, i, 2] = [np.log10(j[5]) * 0.05
                                     for j in fdata[i::5]]
            data_matrix[:, i, 0] = [j[6] + np.random.randn() * 2
                                    for j in fdata[i::5]]
            error_matrix[:, i, 0] = [2 for j in fdata[i::5]]
            data_matrix[:, i, 3] = [np.log10(j[7]) * (1 + np.random.randn() * 0.05)
                                    for j in fdata[i::5]]
            error_matrix[:, i, 3] = [np.log10(j[7]) * 0.05
                                     for j in fdata[i::5]]
            data_matrix[:, i, 1] = [j[8] + 180 + np.random.randn() * 2
                                    for j in fdata[i::5]]
            error_matrix[:, i, 1] = [2 for j in fdata[i::5]]
        for site in fdata[::5]:
            self.x_locs.append(site[0])
            self.y_locs.append(site[1])
        self.data = data_matrix
        self.data_errors = error_matrix

    def evaluate(self, lo_blob_parameters):
        inputs, rms = [], []
        self.lo_blob_parameters = lo_blob_parameters
        for blobs in lo_blob_parameters:
            for idx, (x, y) in enumerate(zip(self.x_locs, self.y_locs)):
                neural_in = (np.array([scale_loc(x), scale_loc(y)] +
                                      blobs.tolist()).reshape(len(blobs) + 2))
                inputs.append(neural_in)
        responses = self.forward(np.array(inputs))
        self.responses = responses
        self.inputs = inputs
        chi_2 = (responses - self.data)**2 / self.data_errors**2
        rms = [np.sqrt(np.mean(chi_2[i]))
               for i in range(len(lo_blob_parameters))]
        return rms

    def forward(self, inputs):
        self.inputs = scale_params(inputs)
        net = self.neural_nets[inputs.shape[1]]
        raw_response = net['model'].predict(self.inputs,
                                            batch_size=self.batch_size)
        scale_response = raw_response * (net['max'] - net['min']) + net['min']
        scale_response[:, 10:] = scale_response[:, 10:]
        scale_response[:, :10] = 90 - scale_response[:, :10]
        # print(self.data.shape)
        print(scale_response[0])
        neural_matrix = np.zeros([len(self.lo_blob_parameters),
                                  self.data.shape[0], self.data.shape[1],
                                  self.data.shape[2]])
        scale_response = scale_response.reshape(len(self.lo_blob_parameters),
                                                -1)
        for bdx, (blobs) in enumerate(self.lo_blob_parameters):
            for idx in range(59):
                neural = scale_response[bdx, idx *
                                        20:(idx + 1) * 20].reshape(2, 5, 2)
                # print neural.shape, neural_matrix.shape
                # print(scale_response[bdx, idx *
                #                         20:(idx + 1) * 20].reshape(2, 5, 2)[0])
                for jdx in range(5):
                    # print(neural_matrix[bdx, idx, jdx].shape)
                    neural_matrix[bdx, idx, jdx] = ([neural[0, jdx, 0],
                                                     neural[0, jdx, 1],
                                                     neural[1, jdx, 0],
                                                     neural[1, jdx, 1]])
        self.raw_response = raw_response
        self.scale_response = neural_matrix
        return neural_matrix


def scale_params(parameters, unscale=False):
    try:
        mins = np.array([0, 0, 0.1, 0, 0.01, 0.02,
                         0.4, 0.4, 0.05, 0.05, 0, 0.05])
        maxs = np.array([1, 1, 0.9, 1, 0.99, 0.99,
                         0.6, 0.6, 0.2, 0.2, 0.4, 0.2])
        return (parameters * (maxs - mins) + mins)
    except:
        mins = mins[2:]
        maxs = maxs[2:]
        return (parameters * (maxs - mins) + mins)


def evol(params, sigma, popsize, maxiter, forward_fn):
    trajectory = []
    cmalogfilenameprefix = 'outcmaes'
    if len(sys.argv) > 1:
        cmalogfilenameprefix = sys.argv[1] + cmalogfilenameprefix
    es = cma.CMAEvolutionStrategy(params, sigma, {
        'popsize': popsize, 'maxiter': maxiter, 'bounds': [0, 1]})
    fit = [256] * popsize
    while not es.stop():
        X = es.ask(popsize)
        fit = forward_fn.evaluate(X)
        es.tell(X, fit)
        es.disp()
        trajectory.append(es.best.get()[:2])
    return trajectory[-1][0], trajectory[-1][1]


def pad_model(b, no_blobs):
    print b
    # no_blobs = len(b[1:])/9
    b = list(b)
    # b[6] = 0
    inserts = []
    for i in range(no_blobs):
        for j in [7, 10, 11]:
            inserts.append(i * 9 + 1 + j)
    for i in inserts[::-1]:
        b.insert(i, 0)
    return b


def get_model(b, model_name='./targetta'):
    print b
    b = pad_model(scale_params(b), 1)
    c = ' '.join(map(str, b))
    print(c)

    model_template = './template'
    target_model = model_name

    callString = './saveModel ' + model_template + ' ' + c + ' > ' + target_model
    os.system(callString)
    time.sleep(2)
    # mtNS = ws2vtk(target_model, '../oneblob.resp')
    # return mtNS

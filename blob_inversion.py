import sys
import keras
import numpy
import os, sys
import numpy as np
import time
from glob import glob1
import cma

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
                 batch_size=1024 * 4, override_error=False, scale=False,
                 read_data=True):
        """
            reads in data file from current directory and neural network data
            it is assumed that neural net models are structured as follows:
                base_name.model -- filename for keras hdf5 model
                base_name.min   -- filename for minimum scale parameters
                base_name.max   -- filename for maximum scale parameters

            currently only programmed for a datafile which exactly matches the
            40 parameter response files which the neural networks have been
        /    trained on
        """

        self.override_error = override_error
        if read_data:
            self.data, self.data_errors = self.read_data(
                data_fn, scale=False, error=error)
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
    def pop_data(self, filename, error=None, off_diag=False):
        with open(filename) as f:
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
            if 'Iteration' in line or 'ERROR' in line:
                break
            if 'DATA' in line:
                period = float(line.split(' ')[-1])
                response[period] = []
                continue
            tensor = [float(i.strip()) for i in line.split(' ') if i.strip()]
            if len(tensor) != 4:
                tensor = tensor[2:6]
            app_xy = np.log10(np.linalg.norm(np.array(tensor[0:2])*796)**2*period*0.2)
            phs_xy = np.degrees(np.arctan2(tensor[0],-tensor[1]))
            app_yx = np.log10(np.linalg.norm(np.array(tensor[2:4])*796)**2*period*0.2)
            phs_yx = np.degrees(np.arctan2(-tensor[2],tensor[3]))
            response[period].append([app_xy, phs_xy, app_yx, phs_yx])
        ins = []
        x_set = sorted(list(set(x_locs)))
        y_set = sorted(list(set(y_locs)))
        periods = sorted(response.keys())
        data_matrix = np.zeros([len(x_set), len(y_set), len(response.keys()), 4])
        for idx, (i, j) in enumerate(zip(x_locs, y_locs)):
            x = np.argmax(np.array(x_set) == i)
            y = np.argmax(np.array(y_set) == j)
            for jdx, k in enumerate(sorted(response.keys())):
                data_matrix[x, y, jdx] = response[k][idx]
        self.x_set = x_set
        self.y_set = y_set
        self.x_locs = x_locs
        self.y_locs = y_locs
        self.periods = periods
        return data_matrix
    def read_data(self, data_fn, scale=False, synth=True, error=0.05):
        data = self.pop_data(data_fn, error=error)
        errors = np.ones(data.shape)
        errors[:, :, :, 0] = errors[:, :, :, 2] = np.arctan(error)*90/np.pi*2
        errors[:, :, :, 1] = data[:, :, :, 1]*error*2
        errors[:, :, :, 3] = data[:, :, :, 3]*error*2
        return data, errors
    def evaluate(self, lo_blob_parameters):
        inputs, rms = [], []
        self.lo_blob_parameters = lo_blob_parameters
        for blobs in lo_blob_parameters:
            for idx, (i, j) in enumerate(zip(self.x_locs, self.y_locs)):
                x = np.argmax(np.array(self.x_set) == i)
                y = np.argmax(np.array(self.y_set) == j)
                neural_in = np.array([scale_loc(self.x_set[x]), scale_loc(self.y_set[y])] + blobs.tolist()).reshape(len(blobs)+2)
                inputs.append(neural_in)
        responses = self.forward(np.array(inputs))
        self.responses = responses
        self.inputs = inputs
        chi_2 = (responses - self.data)**2 / self.data_errors**2
        rms = [np.sqrt(np.mean(chi_2[i])) for i in range(len(lo_blob_parameters))]
        return rms
    def forward(self, inputs):
        self.inputs = scale_params(inputs)
        net = self.neural_nets[inputs.shape[1]]
        raw_response = net['model'].predict(self.inputs, batch_size=self.batch_size)
        scale_response = raw_response * (net['max'] - net['min']) + net['min']
        scale_response[:, 10:] = scale_response[:, 10:]
        scale_response[:, :10] = 90 - scale_response[:, :10]
        neural_matrix = np.zeros([len(self.lo_blob_parameters), self.data.shape[0], self.data.shape[1], self.data.shape[2], self.data.shape[3]])
        scale_response = scale_response.reshape(len(self.lo_blob_parameters), -1)
        for bdx, (blobs) in enumerate(self.lo_blob_parameters):
            for idx, (i, j) in enumerate(zip(self.x_locs, self.y_locs)):
                x = np.argmax(np.array(self.x_set) == i)
                y = np.argmax(np.array(self.y_set) == j)
                neural = scale_response[bdx, idx*20:(idx+1)*20].reshape(2, 5, 2)
                for jdx in range(len(self.periods)):
                    neural_matrix[bdx, x, y, jdx] = ([neural[1, jdx, 0], neural[0, jdx, 0], neural[1, jdx, 1], neural[0, jdx, 1]])
        self.raw_response = raw_response
        self.scale_response = neural_matrix
        return neural_matrix

def scale_params(parameters, unscale=False):
    try:
        mins = np.array([0, 0, 0.1, 0, 0.01, 0.02, 0.4, 0.4, 0.05, 0.05, 0, 0.05])
        maxs = np.array([1, 1, 0.9, 1, 0.99, 0.99, 0.6, 0.6, 0.2, 0.2, 0.4, 0.2])
        return (parameters*(maxs-mins) + mins)
    except:
        mins = mins[2:]
        maxs = maxs[2:]
        return (parameters*(maxs-mins) + mins)

def evol(params, sigma, popsize, maxiter, forward_fn):
    trajectory = []
    es = cma.CMAEvolutionStrategy(params, sigma, {
        'popsize': popsize, 'maxiter': maxiter, 'bounds': [0, 1]})
    fit = [256] * popsize
    while not es.stop():
        X = es.ask(popsize)
        fit = forward_fn.evaluate(X)
        es.tell(X, fit)
        es.disp()
        trajectory.append(es.best.get()[:2])
    return es, trajectory

def brute_guess(params, forward_fn, batch_size, iterations):
    try:
        best_guess, best_x = np.inf, []
        for i in range(iterations):
            guesses = np.random.rand(batch_size, params)
            results = forward_fn.evaluate(guesses)
            if min(results) < best_guess:
                best_guess = min(results)
                best_x = guesses[np.argmin(results)]
            print i, best_guess
        return best_x
    except KeyboardInterrupt:
        return best_x

def swarm(params, forward_fn, swarm_size, phi_p, phi_g, max_iters,
          current_swarm=None):
    start_time = time.time()
    dhigh = 1e-2
    dlow = 5e-4
    if params == '1':
        params = 10
    elif params == '2':
        params = 19
    elif params == '3':
        params = 28
    elif params == '4':
        params = 37
    try:
        current_swarm.shape
    except AttributeError:
        current_swarm = np.random.rand(swarm_size, params)
    bests, i = [], 0
    velocity = np.random.rand(swarm_size, params) * 2 - 1
    best_swarm = current_swarm
    best_results = np.array(forward_fn.evaluate(
        current_swarm)).reshape(swarm_size, 1)
    global_best = best_swarm[np.argmax(best_results)]
    while i < max_iters:
        try:
            i += 1
            omega = 0.3 + 0.4*(max_iters - i + 1)/max_iters
            rp = np.random.rand(params)
            rg = np.random.rand(params)
            diversity = np.min(np.sqrt(np.sum((global_best-current_swarm)**2,axis=0)))
            global_attract = 1 if diversity > dhigh else -1
            local_attract = 1 if diversity > dlow else -1
            velocity = omega * velocity + phi_p * rp * (best_swarm - current_swarm) * local_attract + \
                phi_g * rg * (global_best - current_swarm) * global_attract
            current_swarm = (current_swarm + velocity)
            vel_change = np.random.rand(swarm_size, params)
            velocity = (current_swarm > 0) * velocity + \
                (current_swarm <= 0) * velocity * -vel_change
            velocity = (current_swarm < 1) * velocity + \
                (current_swarm >= 1) * velocity * -vel_change
            current_swarm = current_swarm.clip(0, 1)
            results = np.array(forward_fn.evaluate(
                current_swarm)).reshape(swarm_size, 1)
            new_results = (best_results > results)
            best_results = (best_results < results) * best_results + \
                (best_results >= results) * results
            best_swarm = (best_results < results) * best_swarm + \
                (best_results >= results) * current_swarm
            global_best = best_swarm[np.argmin(best_results)]
            info = 'it: {} RMS {} div {}, {}, {}, {}'.format(i, round(min(best_results)[0],5), diversity, local_attract, global_attract, omega)
            sys.stdout.flush()
            sys.stdout.write("\r" + info)
            bests.append(min(best_results)[0])
        except KeyboardInterrupt:
            elapsed = int(time.time() - start_time)
            info = 'search finished after {}m{}s elapsed with rms of {}'.format(elapsed / 60,
                                                                 elapsed % 60, round(bests[-1],2))
            print(info)
            break
    return {'best_swarm': best_swarm, 'global_best': global_best,
              'current_swarm': current_swarm, 'velocity': velocity, 'bests': bests}

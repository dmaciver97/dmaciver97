import numpy as np
import ultranest
import ultranest.stepsampler
from ultranest.plot import cornerplot

class fit(object):
    def __init__(self, params, x_data, y_data):
        self.param_names = params
        self.x = x_data
        self.y = y_data
            
    def prior_transfrom(self, cube):
        self.params = cube.copy()
        for i in range(len(cube)):
            self.params[i][0] = 10**(cube[i][0]*(8-6)+6)
            self.params[i][1] = 10**(cube[i][1]*(8-6)+6)

        return self.params

    def log_likelihood(self, params):
        likelihood = np.ndarray(shape=(len(params)))

        for i in range(len(params)):
            A = self.params[i][0]
            B = self.params[i][1]
            y_model = [1/(A+B*x**2) for x in self.x]
            likelihood[i] = sum([((y-y_m)/0.05)**2 for y, y_m in zip(self.y, y_model)])

        return likelihood

    def run(self):
        sampler = ultranest.ReactiveNestedSampler(self.param_names,
        self.log_likelihood,
        self.prior_transfrom,
        resume = 'overwrite',
        log_dir = './results',
        vectorized = True)

        print('Start')
        sampler.stepsampler=ultranest.stepsampler.SliceSampler(nsteps=20, 
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction)
        results=sampler.run(Lepsilon = 0.1)
        sampler.print_results()
        sampler.plot_run()
        sampler.plot_corner()
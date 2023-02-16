import numpy
import os
import pandas
import pickle

from os import path
from gpcam.autonomous_experimenter import AutonomousExperimenterGP


class ML_GPCAM_plus_back:

    def __init__(self, fitobj=None, par=None, configuration=None, qmin=0.001, qmax=0.8):
        self.fitobj = fitobj

        self.df_par = par
        self.simpar = pandas.DataFrame(self.df_par.loc['value'])
        self.simpar.reset_index(inplace=True)
        self.simpar.columns = ['par', 'value']

        self.configuration = configuration
        self.qmin = qmin
        self.qmax = qmax

        self.spath = os.getcwd()

        self.my_ae = None
        self.my_ae_back = None

        self.gpcam_init_dataset_size = 20
        self.acq_func = 'shannon_ig_vec'

        self.x = None
        self.y = None
        self.v = None
        self.x_back = None
        self.y_back = None
        self.v_back = None
        self.gpiteration = 0
        self.gpcam_iterations = 10000
        self.gpiteration_back = 0
        self.gpcam_step = 10

        self.gpCAMstream_back = None

        if path.isfile(path.join(self.spath, 'gpCAMstream.pkl')):
            with open(path.join(self.spath, 'gpCAMstream.pkl'), 'rb') as file:
                self.gpCAMstream = pickle.load(file)
                self.gpiteration = len(self.gpCAMstream['position'])
        else:
            self.gpCAMstream = {'position': [], 'value': [], 'variance': []}
            # self.gpCAMstream_back = {'position': [], 'value': [], 'variance': []}

    def gpcam_prediction(self, my_ae=None, position=None):
        res = my_ae.gp_optimizer.posterior_mean(position)
        f = res["f(x)"]
        return f

    def gpcam_SANS_instrument(self, data):
        print("This is the current length of the data received by gpCAM: ", len(data))
        print("Suggested by gpCAM: ", data)
        for entry in data:
            pars = entry['position']
            for i in self.simpar.index:
                self.simpar.at[i, 'value'] = pars[i]
            li_data = self.fitobj.fnSimulateData(liConfigurations=self.configuration, qmin=self.qmin, qmax=self.qmax,
                                                 t_total=None, simpar=self.simpar, save_file=False, verbose=False)
            li_data = li_data[0][1]
            value = li_data['I'].to_numpy()
            variance = li_data['dI'].to_numpy()
            entry["value"] = value
            entry['variance'] = variance
            # entry["cost"]  = [np.array([0,0]),entry["position"],np.sum(entry["position"])]
            self.gpCAMstream['position'].append(entry['position'])
            self.gpCAMstream['value'].append(value)
            self.gpCAMstream['variance'].append(variance)
        return data

    def gpcam_SANS_instrument_back(self, data):
        print("This is the current length of the data received by gpCAM: ", len(data))
        print("Suggested by gpCAM: ", data)
        num = len(data)
        for count, entry in enumerate(data):
            # assign last num results from forward calculation to back
            entry["position"] = self.gpCAMstream['value'][count-num]
            entry["value"] = self.gpCAMstream['position'][count-num]
            entry['variance'] = 0
            # entry["cost"]  = [np.array([0,0]),entry["position"],np.sum(entry["position"])]
        return data

    def gpcam_SANS_training(self, num_sans_points=None, sansmin=None, sansmax=None):
        # Using the gpCAM global optimizer, follows the example from the gpCAM website
        # initialization
        # feel free to try different acquisition functions, e.g. optional_acq_func, "covariance", "shannon_ig"
        # note how costs are defined in for the autonomous experimenter

        parlimits = []
        for par in self.df_par:
            parlimits.append([self.df_par[par]['lowerlimit'], self.df_par[par]['upperlimit']])
        parlimits = numpy.array(parlimits)
        numpars = len(parlimits)

        parlimits_back = []
        for i in range(num_sans_points):
            parlimits_back.append([sansmin, sansmax])
        parlimits_back = numpy.array(parlimits_back)
        numpars_back = num_sans_points

        self.x = self.gpCAMstream['position']
        self.y = self.gpCAMstream['value']
        self.v = self.gpCAMstream['variance']

        if len(self.x) >= 1:
            # use any previously computed results
            self.x = numpy.array(self.x)
            self.y = numpy.array(self.y)
            self.v = numpy.array(self.v)
            self.gpiteration = len(self.x)
            bFirstEval = False
        else:
            self.x = None
            self.y = None
            self.v = None
            self.gpiteration = 0
            bFirstEval = True

        hyperpars = numpy.ones([numpars + 1])
        hyperpars_back = numpy.ones([numpars + 1])
        # the zeroth hyper bound is associated with a signal variance
        # the others with the length scales of the parameter inputs
        hyper_bounds = numpy.array([[0.001, 100]] * (numpars + 1))
        hyper_bounds_back = numpy.array([[0.001, 100]] * (numpars_back + 1))

        for i in range(len(parlimits)):
            delta = parlimits[i][1] - parlimits[i][0]
            hyper_bounds[i + 1] = [delta * 1e-3, delta * 1e1]

        for i in range(len(parlimits_back)):
            delta = parlimits_back[i][1] - parlimits_back[i][0]
            hyper_bounds_back[i + 1] = [delta * 1e-3, delta * 1e1]

        self.my_ae = AutonomousExperimenterGP(parlimits, hyperpars, hyper_bounds,
                                              init_dataset_size=self.gpcam_init_dataset_size,
                                              instrument_func=self.gpcam_SANS_instrument,
                                              acq_func=self.acq_func,  # optional_acq_func,
                                              # cost_func = optional_cost_function,
                                              # cost_update_func = optional_cost_update_function,
                                              x=self.x, y=self.y, v=self.v,
                                              # cost_func_params={"offset": 5.0, "slope": 10.0},
                                              kernel_func=None, use_inv=True,
                                              communicate_full_dataset=False, ram_economy=True)

        self.x_back = self.gpCAMstream['value']
        self.y_back = self.gpCAMstream['position']
        self.v_back = numpy.zeros_like(self.y_back)
        self.gpiteration_back = len(self.x_back)
        self.my_ae_back = AutonomousExperimenterGP(parlimits_back, hyperpars_back, hyper_bounds_back,
                                                   init_dataset_size=0,
                                                   instrument_func=self.gpcam_SANS_instrument_back,
                                                   acq_func=self.acq_func,  # optional_acq_func,
                                                   # cost_func = optional_cost_function,
                                                   # cost_update_func = optional_cost_update_function,
                                                   x=self.x_back, y=self.y_back, v=self.v_back,
                                                   # cost_func_params={"offset": 5.0, "slope": 10.0},
                                                   kernel_func=None, use_inv=True,
                                                   communicate_full_dataset=False, ram_economy=True)

        # save and evaluate initial data set if it has been freshly calculate
        if bFirstEval:
            self.gpcam_save_results()
            self.gpcam_prediction(self.my_ae)

        while len(self.my_ae.x) < self.gpcam_iterations:
            print("length of the dataset: ", len(self.my_ae.x))
            self.my_ae.train(method="global", max_iter=10000)
            self.my_ae_back.train(method="global", max_iter=10000)
            # or not, or both, choose "global","local" and "hgdl"
            # update hyperparameters in case they are optimized asynchronously
            self.my_ae.train(method="local")
            self.my_ae_back.train(method="local")
            # or not, or both, choose between "global","local" and "hgdl"
            # training and client can be killed if desired and in case they are optimized asynchronously
            # self.my_ae.kill_training()
            if self.gpcam_step is not None:
                target_iterations = len(self.my_ae.x) + self.gpcam_step
                retrain_async_at = []
            else:
                target_iterations = self.gpcam_iterations
                retrain_async_at = numpy.logspace(start=numpy.log10(len(self.my_ae.x)),
                                                  stop=numpy.log10(self.gpcam_iterations / 2), num=3, dtype=int)
            # run the autonomous loop
            self.my_ae.go(N=target_iterations,
                          retrain_async_at=retrain_async_at,
                          retrain_globally_at=[],
                          retrain_locally_at=[],
                          acq_func_opt_setting=lambda number: "global" if number % 2 == 0 else "local",
                          training_opt_max_iter=20,
                          training_opt_pop_size=10,
                          training_opt_tol=1e-6,
                          acq_func_opt_max_iter=20,
                          acq_func_opt_pop_size=20,
                          acq_func_opt_tol=1e-6,
                          number_of_suggested_measurements=1,
                          acq_func_opt_tol_adjust=0.1)

            self.my_ae_back.go(N=target_iterations,
                               retrain_async_at=retrain_async_at,
                               retrain_globally_at=[],
                               retrain_locally_at=[],
                               acq_func_opt_setting=lambda number: "global" if number % 2 == 0 else "local",
                               training_opt_max_iter=20,
                               training_opt_pop_size=10,
                               training_opt_tol=1e-6,
                               acq_func_opt_max_iter=20,
                               acq_func_opt_pop_size=20,
                               acq_func_opt_tol=1e-6,
                               number_of_suggested_measurements=1,
                               acq_func_opt_tol_adjust=0.1)

            # training and client can be killed if desired and in case they are optimized asynchronously
            if self.gpcam_step is None:
                self.my_ae.kill_training()
                self.my_ae_back.kill_training()

            self.gpcam_save_results()
            self.gpcam_prediction(self.my_ae)

    def gpcam_save_results(self):
        with open(path.join(self.spath, 'gpCAMstream.pkl'), 'wb') as file:
            pickle.dump(self.gpCAMstream, file)


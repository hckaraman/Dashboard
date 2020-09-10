# region modules
import os
import math
import matplotlib.pyplot as plt
# import nam_fun as nam_f
import numpy as np
# import objectivefunctions as obj
import pandas as pd
import seaborn
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import AnchoredText
from scipy import stats
from scipy.optimize import minimize

# endregion

# pd.plotting.register_matplotlib_converters(explicit=True)
seaborn.set()
np.seterr(all='ignore')

# path = os.path.split(os.path.abspath(__file__))[0]
path = 'C:\\Users\\cagri\\Desktop\\Dashboard\\Test'


def nam_f(x, P, T, E, area, SpinOff):
    c = []
    q = []
    qofmin = 0.4
    beta = 0.1
    pmm = 10
    carea = 1.0
    qofmin = 0.4
    beta = 0.1
    pmm = 10
    deltat = 24
    carea = 1.0
    # csnow = 0
    # Set initial states
    States = np.array([0, 0, 0.9 * x[1], 0, 0, 0, 0, 0.1])
    ## Set initial state conditions
    snow = States[0]
    u = States[1]
    l = States[2]
    if1 = States[3]
    if2 = States[4]
    of1 = States[5]
    of2 = States[6]
    bf = States[7]

    ## Set parameters
    umax = x[0]
    lmax = x[1]
    cqof = x[2]
    ckif = x[3] / deltat
    ck12 = x[4] / deltat
    tof = x[5]
    tif = x[6]
    tg = x[7]
    ckbf = x[8] / deltat
    qs = 0
    csnow = x[9]
    snowtemp = x[10]
    # Extra.area = 440/(3.6*24)
    lfrac = l / lmax
    fact = area
    # SimRR = list()
    Qsim = np.zeros(len(P))
    # f = open("result.txt", "w")
    spin = 0

    for t in range(len(P)):
        # Set boundary conditions
        prec = P[t]
        evap = E[t]
        temp = T[t]
        # --------------- Snow storage and snow melt --------------------
        if temp < snowtemp:
            snow = snow + prec
        else:
            qs = csnow * temp
            if snow < qs:
                qs = snow
                snow = 0
            else:
                snow = snow - qs
        # ---------------- Evapotranspiration module --------------------
        if temp < 0:
            u1 = u
        else:
            u1 = u + prec + qs
        if u1 > evap:
            eau = evap
            eal = 0
        else:
            eau = u1
            eal = (evap - eau) * lfrac

        u2 = min(u1 - eau, umax)

        if (lfrac > tif):
            qif = (lfrac - tif) / (1 - tif) * u2 / ckif
        else:
            qif = 0

        u3 = u1 - eau - qif

        if u3 > umax:
            pn = u3 - umax
            u = umax
        else:
            pn = 0
            u = u3
        # -------------------- Net precipitation ------------------------
        n = int(pn / pmm) + 1
        pnlst = pn - (n - 1) * pmm
        eal = eal / n

        qofsum = 0
        gsum = 0
        # ---------------------------------------------------------------
        for i in range(1, n + 1, 1):
            pn = pmm
            if i == n:
                pn = pnlst
            # ------------------ Overland flow --------------------------

            if lfrac > tof:
                qof = cqof * (lfrac - tof) / (1 - tof) * pn
            else:
                qof = 0

            qofsum = qofsum + qof

            # -------------------- Recharge -----------------------------

            if lfrac > tg:
                g = (lfrac - tg) / (1 - tg) * (pn - qof)
            else:
                g = 0

            gsum = gsum + g

            # --- Lower zone storage ---

            dl = pn - qof - g
            l = l + dl - eal

            if l > lmax:
                gsum = gsum + l - lmax
                l = lmax

            lfrac = l / lmax

        qof = qofsum
        g = gsum
        eal = n * eal
        # ------------------------ Baseflow -----------------------------

        c = math.exp(-1. / ckbf)
        bf = bf * c + g * carea * (1 - c)

        # ------------------------ Interflow ----------------------------

        c = math.exp(-1. / ck12)
        if1 = if1 * c + qif * (1 - c)
        if2 = if2 * c + if1 * (1 - c)
        # ------ Overland flow routing and overland flow component ------

        of = 0.5 * (of1 + of2) / deltat

        if of > qofmin:
            ckqof = ck12 * (of / qofmin) ** (-beta)  ## Elementwise
        else:
            ckqof = ck12

        c = math.exp(-1. / ckqof)

        # ---------------------------------------------------------------

        of1 = of1 * c + qof * (1 - c)
        of2 = of2 * c + of1 * (1 - c)
        # Update state variables
        States[0] = snow
        States[1] = u
        States[2] = l
        States[3] = if1
        States[4] = if2
        States[5] = of1
        States[6] = of2
        States[7] = bf

        # Update simulated value
        if t >= SpinOff:
            Qsim[t] = (fact * (if2 + of2 + bf))
        # s = str(q)
        # f.write(s + '\n')

    return Qsim


def bias(evaluation, simulation):
    """
    Bias as shown in Gupta in Sorooshian (1998), Toward improved calibration of hydrologic models:
    Multiple  and noncommensurable measures of information, Water Resources Research

        .. math::

         Bias=\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Bias
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        bias = np.nansum(obs - sim) / len(obs)
        return float(bias)

    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def pbias(evaluation, simulation):
    """
    Procentual Bias

        .. math::

         PBias= 100 * \\frac{\\sum_{i=1}^{N}(e_{i}-s_{i})}{\\sum_{i=1}^{N}(e_{i})}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: PBias
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        sim = np.array(simulation)
        obs = np.array(evaluation)
        return 100 * (float(np.nansum(sim - obs)) / float(np.nansum(obs)))

    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def nashsutcliffe(evaluation, simulation):
    """
    Nash-Sutcliffe model efficinecy

        .. math::

         NSE = 1-\\frac{\\sum_{i=1}^{N}(e_{i}-s_{i})^2}{\\sum_{i=1}^{N}(e_{i}-\\bar{e})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Nash-Sutcliff model efficiency
    :rtype: float

    """
    if len(evaluation) == len(simulation):
        s, e = np.array(simulation), np.array(evaluation)
        # s,e=simulation,evaluation
        mean_observed = np.nanmean(e)
        # compute numerator and denominator
        numerator = np.nansum((e - s) ** 2)
        denominator = np.nansum((e - mean_observed) ** 2)
        # compute coefficient
        return 1 - (numerator / denominator)

    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def lognashsutcliffe(evaluation, simulation, epsilon=0):
    """
    log Nash-Sutcliffe model efficiency

        .. math::

         NSE = 1-\\frac{\\sum_{i=1}^{N}(log(e_{i})-log(s_{i}))^2}{\\sum_{i=1}^{N}(log(e_{i})-log(\\bar{e})^2}-1)*-1

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :epsilon: Value which is added to simulation and evaluation data to errors when simulation or evaluation data has zero values
    :type: float or list

    :return: log Nash-Sutcliffe model efficiency
    :rtype: float

    """
    if len(evaluation) == len(simulation):
        s, e = np.array(simulation) + epsilon, np.array(evaluation) + epsilon
        return float(1 - sum((np.log(s) - np.log(e)) ** 2) / sum((np.log(e) - np.mean(np.log(e))) ** 2))
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def log_p(evaluation, simulation):
    """
    Logarithmic probability distribution

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Logarithmic probability distribution
    :rtype: float
    """
    scale = np.mean(evaluation) / 10
    if scale < .01:
        scale = .01
    if len(evaluation) == len(simulation):
        y = (np.array(evaluation) - np.array(simulation)) / scale
        normpdf = -y ** 2 / 2 - np.log(np.sqrt(2 * np.pi))
        return np.mean(normpdf)
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def correlationcoefficient(evaluation, simulation):
    """
    Correlation Coefficient

        .. math::

         r = \\frac{\\sum ^n _{i=1}(e_i - \\bar{e})(s_i - \\bar{s})}{\\sqrt{\\sum ^n _{i=1}(e_i - \\bar{e})^2} \\sqrt{\\sum ^n _{i=1}(s_i - \\bar{s})^2}}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Corelation Coefficient
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        correlation_coefficient = np.corrcoef(evaluation, simulation)[0, 1]
        return correlation_coefficient
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def rsquared(evaluation, simulation):
    """
    Coefficient of Determination

        .. math::

         r^2=(\\frac{\\sum ^n _{i=1}(e_i - \\bar{e})(s_i - \\bar{s})}{\\sqrt{\\sum ^n _{i=1}(e_i - \\bar{e})^2} \\sqrt{\\sum ^n _{i=1}(s_i - \\bar{s})^2}})^2

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Coefficient of Determination
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        return correlationcoefficient(evaluation, simulation) ** 2
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def mse(evaluation, simulation):
    """
    Mean Squared Error

        .. math::

         MSE=\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Mean Squared Error
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mse = np.nanmean((obs - sim) ** 2)
        return mse
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def rmse(evaluation, simulation):
    """
    Root Mean Squared Error

        .. math::

         RMSE=\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Root Mean Squared Error
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        return np.sqrt(mse(evaluation, simulation))
    else:
        # logging.warning("evaluation and simulation lists do not have the same length.")
        return np.nan


def mae(evaluation, simulation):
    """
    Mean Absolute Error

        .. math::

         MAE=\\frac{1}{N}\\sum_{i=1}^{N}(\\left |  e_{i}-s_{i} \\right |)

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Mean Absolute Error
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        obs, sim = np.array(evaluation), np.array(simulation)
        mae = np.mean(np.abs(sim - obs))
        return mae
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def rrmse(evaluation, simulation):
    """
    Relative Root Mean Squared Error

        .. math::

         RRMSE=\\frac{\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}}{\\bar{e}}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Relative Root Mean Squared Error
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        rrmse = rmse(evaluation, simulation) / np.mean(evaluation)
        return rrmse
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def agreementindex(evaluation, simulation):
    """
    Agreement Index (d) developed by Willmott (1981)

        .. math::

         d = 1 - \\frac{\\sum_{i=1}^{N}(e_{i} - s_{i})^2}{\\sum_{i=1}^{N}(\\left | s_{i} - \\bar{e} \\right | + \\left | e_{i} - \\bar{e} \\right |)^2}


    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Agreement Index
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        simulation, evaluation = np.array(simulation), np.array(evaluation)
        Agreement_index = 1 - (np.sum((evaluation - simulation) ** 2)) / (np.sum(
            (np.abs(simulation - np.mean(evaluation)) + np.abs(evaluation - np.mean(evaluation))) ** 2))
        return Agreement_index
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def covariance(evaluation, simulation):
    """ Covariance

        .. math::
         Covariance = \\frac{1}{N} \\sum_{i=1}^{N}((e_{i} - \\bar{e}) * (s_{i} - \\bar{s}))

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Covariance
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        obs_mean = np.mean(obs)
        sim_mean = np.mean(sim)
        covariance = np.mean((obs - obs_mean) * (sim - sim_mean))
        return covariance
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def decomposed_mse(evaluation, simulation):
    """
    Decomposed MSE developed by Kobayashi and Salam (2000)

        .. math ::
         dMSE = (\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i}))^2 + SDSD + LCS

         SDSD = (\\sigma(e) - \\sigma(s))^2

         LCS = 2 \\sigma(e) \\sigma(s) * (1 - \\frac{\\sum ^n _{i=1}(e_i - \\bar{e})(s_i - \\bar{s})}{\\sqrt{\\sum ^n _{i=1}(e_i - \\bar{e})^2} \\sqrt{\\sum ^n _{i=1}(s_i - \\bar{s})^2}})

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Decomposed MSE
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        e_std = np.std(evaluation)
        s_std = np.std(simulation)

        bias_squared = bias(evaluation, simulation) ** 2
        sdsd = (e_std - s_std) ** 2
        lcs = 2 * e_std * s_std * (1 - correlationcoefficient(evaluation, simulation))

        decomposed_mse = bias_squared + sdsd + lcs

        return decomposed_mse
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def kge(evaluation, simulation, return_all=False):
    """
    Kling-Gupta Efficiency

    Corresponding paper:
    Gupta, Kling, Yilmaz, Martinez, 2009, Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling

    output:
        kge: Kling-Gupta Efficiency
    optional_output:
        cc: correlation
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    if len(evaluation) == len(simulation):
        cc = np.corrcoef(evaluation, simulation)[0, 1]
        alpha = np.std(simulation) / np.std(evaluation)
        beta = np.sum(simulation) / np.sum(evaluation)
        kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        if return_all:
            return kge, cc, alpha, beta
        else:
            return kge
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def _spearmann_corr(x, y):
    """Separmann correlation coefficient"""
    col = [list(a) for a in zip(x, y)]
    xy = sorted(col, key=lambda x: x[0], reverse=False)
    # rang of x-value
    for i, row in enumerate(xy):
        row.append(i + 1)

    a = sorted(xy, key=lambda x: x[1], reverse=False)
    # rang of y-value
    for i, row in enumerate(a):
        row.append(i + 1)

    MW_rank_x = np.nanmean(np.array(a)[:, 2])
    MW_rank_y = np.nanmean(np.array(a)[:, 3])

    numerator = np.nansum([float((a[j][2] - MW_rank_x) * (a[j][3] - MW_rank_y)) for j in range(len(a))])
    denominator1 = np.sqrt(np.nansum([(a[j][2] - MW_rank_x) ** 2. for j in range(len(a))]))
    denominator2 = np.sqrt(np.nansum([(a[j][3] - MW_rank_x) ** 2. for j in range(len(a))]))
    return float(numerator / (denominator1 * denominator2))


def kge_non_parametric(evaluation, simulation, return_all=False):
    """
    Non parametric Kling-Gupta Efficiency

    Corresponding paper:
    Pool, Vis, and Seibert, 2018 Evaluating model performance: towards a non-parametric variant of the Kling-Gupta efficiency, Hydrological Sciences Journal.

    output:
        kge: Kling-Gupta Efficiency

    author: Nadine Maier and Tobias Houska
    optional_output:
        cc: correlation
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    if len(evaluation) == len(simulation):
        ## self-made formula
        cc = _spearmann_corr(evaluation, simulation)

        ### scipy-Version
        # cc = stm.spearmanr(evaluation, simulation, axis=0)[0]

        ### pandas version
        # a  = pd.DataFrame({'eva': evaluation, 'sim': simulation})
        # cc = a.ix[:,1].corr(a.ix[:,0], method = 'spearman')

        fdc_sim = np.sort(simulation / (np.nanmean(simulation) * len(simulation)))
        fdc_obs = np.sort(evaluation / (np.nanmean(evaluation) * len(evaluation)))
        alpha = 1 - 0.5 * np.nanmean(np.abs(fdc_sim - fdc_obs))

        beta = np.mean(simulation) / np.mean(evaluation)
        kge = 1 - np.sqrt((cc - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        if return_all:
            return kge, cc, alpha, beta
        else:
            return kge
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def rsr(evaluation, simulation):
    """
    RMSE-observations standard deviation ratio

    Corresponding paper:
    Moriasi, Arnold, Van Liew, Bingner, Harmel, Veith, 2007, Model Evaluation Guidelines for Systematic Quantification of Accuracy in Watershed Simulations

    output:
        rsr: RMSE-observations standard deviation ratio
    """
    if len(evaluation) == len(simulation):
        rsr = rmse(evaluation, simulation) / np.std(evaluation)
        return rsr
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


def volume_error(evaluation, simulation):
    """
    Returns the Volume Error (Ve).
    It is an indicator of the agreement between the averages of the simulated
    and observed runoff (i.e. long-term water balance).
    used in this paper:
    Reynolds, J.E., S. Halldin, C.Y. Xu, J. Seibert, and A. Kauffeldt. 2017.
    “Sub-Daily Runoff Predictions Using Parameters Calibrated on the Basis of Data with a
    Daily Temporal Resolution.” Journal of Hydrology 550 (July):399–411.
    https://doi.org/10.1016/j.jhydrol.2017.05.012.

        .. math::

         Sum(simulation-evaluation)/sum(simulation)
    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Volume Error
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        ve = np.sum(simulation - evaluation) / np.sum(evaluation)
        return float(ve)
    else:
        # logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


_all_functions = [agreementindex, bias, correlationcoefficient, covariance, decomposed_mse,
                  kge, log_p, lognashsutcliffe, mae, mse, nashsutcliffe, pbias, rmse, rrmse, rsquared,
                  rsr, volume_error
                  ]


def calculate_all_functions(evaluation, simulation):
    """
    Calculates all objective functions from spotpy.objectivefunctions
    and returns the results as a list of name/value pairs

    :param evaluation: a sequence of evaluation data
    :param simulation: a sequence of simulation data
    :return: A list of (name, value) tuples
    """

    result = []
    for f in _all_functions:
        # Check if the name is not private and attr is a function but not this

        try:
            result.append((f.__name__, f(evaluation, simulation)))
            text = f.__name__ + " = "
            # print(f.__name__ + " =  %f" %f(evaluation, simulation) )
            # print("Root Mean Square Error (RMSE) = %f" % RMSE)

        except:
            result.append((f.__name__, np.nan))

    return result


class Nam(object):
    _dir = r'D:\DRIVE\TUBITAK\Hydro_Model\Data\Darbogaz'
    _data = "Darbogaz.csv"

    def __init__(self, area, input_parameters, calibration=False, objective='NSE'):
        self._working_directory = None
        self.Data_file = None
        self.df = None
        self.P = None
        self.T = None
        self.E = None
        self.Qobs = None
        self.area = area / (3.6 * 24)
        self.Area = area
        self.Spinoff = 0
        self.parameters = None
        self.Qfit = None
        self.dfh = None
        self.initial_cal = np.array([10, 100, 0.5, 500, 10, 0.5, 0.5, 0, 2000, 2.15, 2])
        # self.initial = np.array([5.59441567e+00,6.85168038e+02,1.30412167e-01,8.47239393e+02,4.00934557e+01,4.21557738e-01,4.88201564e-01,4.09627612e-02,1.67517734e+03,4.09537018e-01,3.71693424e+00])
        self.initial = np.array(input_parameters)
        self.Qsim = None
        self.n = None
        self.Date = None
        self.bounds = (
            (0.01, 50), (0.01, 1000), (0.01, 1), (200, 1000), (10,
                                                               50), (0.01, 0.99), (0.01, 0.99), (0.01, 0.99),
            (500, 5000), (0, 4), (-2, 4))
        self.NSE = None
        self.RMSE = None
        self.PBIAS = None
        self.Cal = calibration
        self.statistics = None
        self.export = 'Result.csv'
        self.flowduration = None
        self.objective = objective

    @property
    def process_path(self):
        return self._working_directory

    @process_path.setter
    def process_path(self, value):
        self._working_directory = value
        pass

    def DataRead(self):
        self.df = pd.read_csv(self.Data_file, sep=',',
                              parse_dates=[0], header=0)
        self.df = self.df.set_index('Date')

    def InitData(self):
        self.P = self.df.P
        self.T = self.df.Temp
        self.E = self.df.E
        self.Qobs = self.df.Q
        self.n = self.df.__len__()
        self.Qsim = np.zeros(self.n)
        self.Date = self.df.index.to_pydatetime()

    def nash(self, qobserved, qsimulated):
        s, e = np.array(qobserved), np.array(qsimulated)
        # s,e=simulation,evaluation
        mean_observed = np.nanmean(e)
        # compute numerator and denominator
        numerator = np.nansum((e - s) ** 2)
        denominator = np.nansum((e - mean_observed) ** 2)
        # compute coefficient
        return 1 - (numerator / denominator)

    def Objective(self, x):
        self.Qsim = nam_f(
            x, self.P, self.T, self.E, self.area, self.Spinoff)
        # n = math.sqrt((sum((self.Qsim - self.Qobs) ** 2)) / len(self.Qobs))
        # n = self.nash(self.Qobs, self.Qsim)
        # mean_observed = np.nanmean(self.Qobs)
        # numerator = np.nansum((self.Qobs - self.Qsim) ** 2)
        # denominator = np.nansum((self.Qobs - mean_observed) ** 2)
        # n = 1 - (numerator / denominator)
        # Q = np.where(self.Qobs > 10, np.nan, self.Qobs)
        # n = obj.nashsutcliffe(self.Qobs,self.Qsim)
        if self.objective == 'RMSE':
            n = rmse(self.Qobs, self.Qsim)
        elif self.objective == 'KGE':
            n = kge(self.Qobs, self.Qsim)
            n = 1 - n
        elif self.objective == 'MAE':
            n = mae(self.Qobs, self.Qsim)
        elif self.objective == 'Volume Error':
            n = volume_error(self.Qobs, self.Qsim)
        else:
            n = nashsutcliffe(self.Qobs, self.Qsim)
            n = 1 - n

        return n

    def run(self):
        self.DataRead()
        self.InitData()
        if self.Cal:
            self.parameters = minimize(self.Objective, self.initial_cal, method='SLSQP', bounds=self.bounds,
                                       options={'maxiter': 1e8, 'disp': False})
            self.Qsim = nam_f(
                self.parameters.x, self.P, self.T, self.E, self.area, self.Spinoff)
            self.parameters = self.parameters.x

        else:
            self.Qsim = nam_f(
                self.initial, self.P, self.T, self.E, self.area, self.Spinoff)
            self.parameters = self.initial
        self.df['Qsim'] = self.Qsim
        return self.df, self.parameters

    def update(self):
        # fit = self.interpolation()
        # self.Qfit = fit(self.Qobs)
        self.df['Qsim'] = self.Qsim
        # self.df['Qfit'] = self.Qfit
        self.flowduration = pd.DataFrame()
        self.flowduration['Qsim_x'] = self.flowdur(self.Qsim)[0]
        self.flowduration['Qsim_y'] = self.flowdur(self.Qsim)[1]
        self.flowduration['Qobs_x'] = self.flowdur(self.Qobs)[0]
        self.flowduration['Qobs_y'] = self.flowdur(self.Qobs)[1]
        # self.df.to_csv(os.path.join(self.process_path, self.export), index=True, header=True)

    def stats(self):
        mean = np.mean(self.Qobs)
        # mean2 = np.mean(self.Qsim)
        self.NSE = 1 - (np.nansum((self.Qsim - self.Qobs) ** 2) /
                        np.nansum((self.Qobs - mean) ** 2))
        self.RMSE = np.sqrt(np.nansum((self.Qsim - self.Qobs) ** 2) / len(self.Qsim))
        self.PBIAS = (np.nansum(self.Qobs - self.Qsim) / np.nansum(self.Qobs)) * 100
        self.statistics = calculate_all_functions(self.Qobs, self.Qsim)

    def interpolation(self):

        idx = np.isfinite(self.Qobs) & np.isfinite(self.Qsim)
        fit = np.polyfit(self.Qobs[idx], self.Qsim[idx], 1)
        fit_fn = np.poly1d(fit)
        return fit_fn

    def draw(self):
        self.stats()
        fit = self.interpolation()
        self.Qfit = fit(self.Qobs)
        width = 15  # Figure width
        height = 10  # Figure height
        f = plt.figure(figsize=(width, height))
        widths = [2, 2, 2]
        heights = [2, 3, 1]
        gs = GridSpec(3, 3, figure=f, width_ratios=widths,
                      height_ratios=heights)
        ax1 = f.add_subplot(gs[1, :])
        ax2 = f.add_subplot(gs[0, :], sharex=ax1)
        ax3 = f.add_subplot(gs[-1, 0])
        ax4 = f.add_subplot(gs[-1, -1])
        ax5 = f.add_subplot(gs[-1, -2])
        color = 'tab:blue'
        ax2.set_ylabel('Precipitation ,mm ', color=color,
                       style='italic', fontweight='bold', labelpad=20, fontsize=13)
        ax2.bar(self.Date, self.df.P, color=color,
                align='center', alpha=0.6, width=1)
        ax2.tick_params(axis='y', labelcolor=color)
        # ax2.set_ylim(0, max(self.df.P) * 1.1, )
        ax2.set_ylim(max(self.df.P) * 1.1, 0)
        ax2.legend(['Precipitation'])
        color = 'tab:red'
        ax2.set_title('NAM Simulation', style='italic',
                      fontweight='bold', fontsize=16)
        ax1.set_ylabel(r'Discharge m$^3$/s', color=color,
                       style='italic', fontweight='bold', labelpad=20, fontsize=13)
        ax1.plot(self.Date, self.Qobs, 'b-', self.Date,
                 self.Qsim, 'r--', linewidth=2.0)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.tick_params(axis='x', labelrotation=45)
        ax1.set_xlabel('Date', style='italic',
                       fontweight='bold', labelpad=20, fontsize=13)
        ax1.legend(('Observed Run-off', 'Simulated Run-off'), loc=2)
        plt.setp(ax2.get_xticklabels(), visible=False)
        anchored_text = AnchoredText("NSE = %.2f\nRMSE = %0.2f\nPBIAS = %0.2f" % (self.NSE, self.RMSE, self.PBIAS),
                                     loc=1, prop=dict(size=11))
        ax1.add_artist(anchored_text)
        # plt.subplots_adjust(hspace=0.05)
        ax3.set_title('Flow Duration Curve', fontsize=11, style='italic')
        ax3.set_yscale("log")
        ax3.set_ylabel(r'Discharge m$^3$/s', style='italic',
                       fontweight='bold', labelpad=20, fontsize=9)
        ax3.set_xlabel('Percentage Exceedence (%)', style='italic',
                       fontweight='bold', labelpad=20, fontsize=9)
        exceedence, sort, low_percentile, high_percentile = self.flowdur(
            self.Qsim)
        ax3.legend(['Precipitation'])
        ax3.plot(self.flowdur(self.Qsim)[0], self.flowdur(self.Qsim)[1], 'b-', self.flowdur(self.Qobs)[0],
                 self.flowdur(self.Qobs)[1], 'r--')
        # ax3.plot(self.flowdur(self.Qobs)[0], self.flowdur(self.Qobs)[1])
        ax3.legend(('Observed', 'Simulated'),
                   loc="upper right", prop=dict(size=7))

        plt.grid(True, which="minor", ls="-")

        st = stats.linregress(self.Qobs, self.Qsim)
        # ax4.set_yscale("log")
        # ax4.set_xscale("log")
        ax4.set_title('Regression Analysis', fontsize=11, style='italic')
        ax4.set_ylabel(r'Simulated', style='italic',
                       fontweight='bold', labelpad=20, fontsize=9)
        ax4.set_xlabel('Observed', style='italic',
                       fontweight='bold', labelpad=20, fontsize=9)
        anchored_text = AnchoredText("y = %.2f\n$R^2$ = %0.2f" % (
            st[0], (st[2]) ** 2), loc=4, prop=dict(size=7))
        # ax4.plot(self.Qobs, fit(self.Qsim), '--k')
        # ax4.scatter(self.Qsim, self.Qobs)
        ax4.plot(self.Qobs, self.Qsim, 'bo', self.Qobs, self.Qfit, '--k')
        ax4.add_artist(anchored_text)

        self.update()
        self.dfh = self.df.resample('M').mean()
        Date = self.dfh.index.to_pydatetime()
        ax5.set_title('Monthly Mean', fontsize=11, style='italic')
        ax5.set_ylabel(r'Discharge m$^3$/s', color=color,
                       style='italic', fontweight='bold', labelpad=20, fontsize=9)
        # ax5.set_xlabel('Date', style='italic', fontweight='bold', labelpad=20, fontsize=9)
        ax5.tick_params(axis='y', labelcolor=color)
        ax5.tick_params(axis='x', labelrotation=45)
        # ax5.set_xlabel('Date', style='italic', fontweight='bold', labelpad=20, fontsize=9)
        ax5.legend(('Observed', 'Simulated'), loc="upper right")
        exceedence, sort, low_percentile, high_percentile = self.flowdur(
            self.Qsim)
        ax5.tick_params(axis='x', labelsize=9)
        ax5.plot(Date, self.dfh.Q, 'b-', Date,
                 self.dfh.Qsim, 'r--', linewidth=2.0)
        ax5.legend(('Observed', 'Simulated'), prop={'size': 7}, loc=1)
        # ax5.plot(dfh.Q)
        # ax5.plot(dfh.Qsim)
        # ax5.legend()
        plt.grid(True, which="minor", ls="-")
        plt.subplots_adjust(hspace=0.03)
        f.tight_layout()
        plt.show()

    def flowdur(self, x):
        exceedence = np.arange(1., len(np.array(x)) + 1) / len(np.array(x))
        exceedence *= 100
        sort = np.sort(x, axis=0)[::-1]
        low_percentile = np.percentile(sort, 5, axis=0)
        high_percentile = np.percentile(sort, 95, axis=0)
        return exceedence, sort, low_percentile, high_percentile

    def drawflow(self):
        f = plt.figure(figsize=(15, 10))
        ax = f.add_subplot(111)
        # fig, ax = plt.subplots(1, 1)
        ax.set_yscale("log")
        ax.set_ylabel(r'Discharge m$^3$/s', style='italic',
                      fontweight='bold', labelpad=20, fontsize=13)
        ax.set_xlabel('Percentage Exceedence (%)', style='italic',
                      fontweight='bold', labelpad=20, fontsize=13)
        exceedence, sort, low_percentile, high_percentile = self.flowdur(
            self.Qsim)
        ax.plot(self.flowdur(self.Qsim)[0], self.flowdur(self.Qsim)[1])
        ax.plot(self.flowdur(self.Qobs)[0], self.flowdur(self.Qobs)[1])
        plt.grid(True, which="minor", ls="-")
        # ax.fill_between(exceedence, low_percentile, high_percentile)
        # plt.show()
        return ax


# Sample Run

def run(area, params, cal, file,obj):
    # if __name__ == '__main__':
    #     params = [6.96780205e+00, 4.86098809e+02, 6.66247792e-01, 5.42601108e+02
    #         , 2.43815545e+01, 8.21285865e-01, 1.00000000e-02, 1.00000000e-02
    #         , 7.37979357e+02, 9.64180895e-01, 2.06295770e+00]
    # print(area,params,cal)
    # area = 97
    # cal = True
    n = Nam(area, params, calibration=cal, objective=obj)
    n.process_path = path
    n.Data_file = os.path.join(n.process_path, 'Data', file + ".csv")
    n.run()
    n.stats()
    n.update()
    # n.draw()
    return n.df, n.parameters, n.statistics, n.flowduration
# run()

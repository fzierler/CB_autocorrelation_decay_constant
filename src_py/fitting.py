from functools import partial

import corrfitter as cf
import gvar as gv
import logging
import numpy as np
from scipy.optimize import curve_fit, minimize
import warnings

from bootstrap import BootstrapSampleSet

warnings.filterwarnings("ignore")


def make_models(tmin, tmax, tp):
    """Create corrfitter model for G(t)."""
    return [cf.Corr2(datatag="Gab", tp=tp, tmin=tmin, tmax=tmax, a="a", b="a", dE="dE")]


def sigle_state_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior["log(a)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(dE)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    return prior


def first_fit_parameters(fit):
    parameters = fit.p
    Energy = np.cumsum(parameters["dE"])
    a_matrix_element = parameters["a"]
    chi2 = fit.chi2
    dof = fit.dof
    return Energy, a_matrix_element, chi2, dof


def fit_correlator_without_bootstrap(
    data_corr, t_lattice, tmin, tmax, Nmax, tp, p0, plotting=False, printing=False
):
    t_lattice = abs(t_lattice)

    fitter = cf.CorrFitter(models=make_models(tmin, tmax, tp))

    for N in range(1, Nmax + 1):
        prior = sigle_state_prior(N)
        fit = fitter.lsqfit(data=data_corr, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print("nterm =", N, 30 * "=")
            print(fit)

    E, a, chi2, dof = first_fit_parameters(fit)
    if plotting:
        fit.show_plots(view="log")

    return E, a, chi2, dof

def fit_cosh_bootstrap(C, plateau_start, plateau_end):
    """This function fits the correlators with a cosh function"""

    C_boot = C.samples

    num_sample = C_boot.shape[0]
    lattice_t = C_boot.shape[1]

    def func(t, a, M):
        return a * a * M * (np.exp(-M * t) + np.exp(-M * (lattice_t - t))) / 2

    x0, pcov = curve_fit(
        func,
        np.arange(plateau_start, plateau_end),
        C.mean[0, plateau_start:plateau_end],
    )

    p0 = dict(
        {"log(a)": np.array([np.log(abs(x0[0]))]), "log(dE)": np.array([np.log(x0[1])])}
    )

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov = np.cov(C_boot.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(C_boot[n], cov))

        E, a, chi2, dof = fit_correlator_without_bootstrap(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            lattice_t,
            p0,
            plotting=False,
            printing=False,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        chi2_dof[n] = chi2 / dof

    correlator_set = dict(Gab=gv.gvar(C.mean[0], cov))
    E_mean, a_mean, chi2, dof = fit_correlator_without_bootstrap(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        lattice_t,
        p0,
        plotting=False,
        printing=False,
    )

    E_fit = BootstrapSampleSet(gv.mean(E_mean[0]), E_sample)
    A_fit = BootstrapSampleSet(
        gv.mean(a_mean[0]) / np.sqrt(gv.mean(E_mean[0])), a_sample / np.sqrt(E_sample)
    )

    return E_fit, A_fit, chi2 / dof


def simultaneous_model(T, tmin, tmax, tp):
    """Create corrfitter model for Gaa(t) and Gab(t)."""
    model = [
        cf.Corr2(datatag="Gaa", tp=tp, tmin=tmin, tmax=tmax, a="a", b="a", dE="dE"),
        cf.Corr2(datatag="Gab", tp=-tp, tmin=tmin, tmax=tmax, a="a", b="b", dE="dE"),
    ]
    return model


def simultaneous_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior["log(a)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(b)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(dE)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    return prior


def simultaneous_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p["dE"])
    a = p["a"]
    b = p["b"]
    chi2 = fit.chi2
    dof = fit.dof
    return E, a, b, chi2, dof


def fit_correlator_simultaneous(
    data_corrs, T, tmin, tmax, Nmax, tp, p0, plotting=False, printing=False
):
    T = abs(T)

    fitter = cf.CorrFitter(models=simultaneous_model(T, tmin, tmax, tp))
    for N in range(1, Nmax + 1):
        prior = simultaneous_prior(N)
        fit = fitter.lsqfit(data=data_corrs, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print("nterm =", N, 30 * "=")
            print(fit)

    E, a, b, chi2, dof = simultaneous_fit_parameters(fit)

    if plotting:
        fit.show_plots(view="log")

    return E, a, b, chi2, dof


def fit_coshsinh_simultaneous(Corr_ss, Corr_sp, plateau_start, plateau_end, lattice_t):
    """This function fits the correlators with cosh and sinh functions simultaneously"""

    x0 = sim_coshsinh_fit(
        Corr_sp.mean[0], Corr_ss.mean[0], lattice_t, plateau_start, plateau_end
    )

    p0 = dict(
        {
            "log(a)": np.array([np.log(abs(x0[0]))]),
            "log(b)": np.array([np.log(x0[1])]),
            "log(dE)": np.array([np.log(x0[2])]),
        }
    )

    # load the ensamble info
    Css = Corr_ss.samples
    Csp = Corr_sp.samples

    num_sample = Css.shape[0]

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    b_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov_ss = np.cov(Css.T)
    cov_sp = np.cov(Csp.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(Csp[n], cov_sp), Gaa=gv.gvar(Css[n], cov_ss))

        E, a, b, chi2, dof = fit_correlator_simultaneous(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            lattice_t,
            p0,
            plotting=False,
            printing=False,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        b_sample[n] = gv.mean(b[0])
        chi2_dof[n] = chi2 / dof

    correlator_set = dict(
        Gab=gv.gvar(Corr_sp.mean[0], cov_sp), Gaa=gv.gvar(Corr_ss.mean[0], cov_ss)
    )

    E_mean, a_mean, b_mean, chi2, dof = fit_correlator_simultaneous(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        lattice_t,
        p0,
        plotting=False,
        printing=False,
    )

    E_fit = BootstrapSampleSet(gv.mean(E_mean[0]), E_sample)
    A_fit = BootstrapSampleSet(
        gv.mean(b_mean[0]) / np.sqrt(gv.mean(E_mean[0])), b_sample / np.sqrt(E_sample)
    )

    return E_fit, A_fit, chi2 / dof


def sim_coshsinh_fit(C1, C2, T, ti, tf):
    y1 = C1[ti:tf]
    y2 = C2[ti:tf]
    comboY = np.append(y1, y2)
    t_slice = np.arange(ti, tf)
    comboX = np.append(t_slice, t_slice)

    def func_sp(t, As, f, M):
        return As * f * (np.exp(-M * t) - np.exp(-M * (T - t))) / 2.0

    def func_ss(t, As, f, M):
        return As**2 * (np.exp(-M * t) + np.exp(-M * (T - t))) / (2.0 * M)

    def comboFunc(
        comboData, a, b, c
    ):  # single data set passed in, extract separate data
        extract1 = comboData[: len(y1)]  # first data
        extract2 = comboData[len(y2) :]  # second data

        result1 = func_sp(extract1, a, b, c)
        result2 = func_ss(extract2, a, b, c)

        return np.append(result1, result2)

    # curve fit the combined data to the combined function
    fittedParameters, pcov = curve_fit(comboFunc, comboX, comboY)

    return fittedParameters